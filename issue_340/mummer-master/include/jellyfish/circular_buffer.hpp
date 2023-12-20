/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __JELLYFISH_CIRCULAR_BUFFER_HPP__
#define __JELLYFISH_CIRCULAR_BUFFER_HPP__

#include <limits>
#include <atomic>


#include <jellyfish/divisor.hpp>

namespace jflib {
  template<typename T, unsigned int n, T g = ((T)1 << n) - 1>
  class basic_circular_buffer {
    static const unsigned int m = std::numeric_limits<T>::digits - n;
    struct splitT {
      T id:m;
      T val:n;
    };
    union elt {
      T      binary;
      splitT split;
    };

    divisor64            _size;
    std::atomic<T>*      _buffer;
    std::atomic<size_t>  _head;
    std::atomic<size_t>  _tail;
    std::atomic<bool>    _closed;

  public:
    typedef T value_type;
    static const T guard = g;

    basic_circular_buffer(size_t size)
      : _size(size+1)
      , _buffer(new std::atomic<T>[_size.d()])
      , _head(0)
      , _tail(0)
      , _closed(false)
    {
      elt init;
      init.split.id  = 0;
      init.split.val = guard;
      for(size_t i = 0; i < _size.d(); ++i)
        _buffer[i] = init.binary;
    }
    virtual ~basic_circular_buffer() {
      delete [] _buffer;
    }

    /** Enqueue an element.
     * @return false if the FIFO is full.
     */
    bool enqueue(const T &v);
    /** Enqueue an element, optimization. No check is made that the
     * FIFO is full. Undetermined behavior if an element is inserted
     * in a full FIFO.
     */
    void enqueue_no_check(const T &v);
    /** Dequeue an element.
     * @return 0 if the FIFO is empty.
     */
    T dequeue();
    bool is_closed() const { return std::atomic_load(&_closed); }
    void close() { _closed.store(true); }

    /// Return capacity of circular buffer
    size_t size() { return _size.d(); }
    /// Return the number of element currently in circular buffer
    size_t fill() {
      size_t head, tail;
      size_t nhead = _head.load(std::memory_order_acquire);
      do {
        head = nhead;
        tail = _tail.load(std::memory_order_acquire);
      } while(head != (nhead = _head.load(std::memory_order_acquire)));

      return head >= tail ? head - tail : head + _size.d() - tail;
    }
  };

  template<typename T, T g = (T)-1>
  class circular_buffer : public basic_circular_buffer<uint64_t, std::numeric_limits<T>::digits, g> {
  public:
    circular_buffer(size_t size) :
      basic_circular_buffer<uint64_t, std::numeric_limits<T>::digits, g>(size) { }
    virtual ~circular_buffer() { }

    bool enqueue(const T &v) {
      return basic_circular_buffer<uint64_t, std::numeric_limits<T>::digits, g>::enqueue((uint64_t)v);
    }
    T dequeue() {
      return basic_circular_buffer<uint64_t, std::numeric_limits<T>::digits, g>::dequeue();
    }
  };
}

template<typename T, unsigned int n, T guard>
bool jflib::basic_circular_buffer<T,n,guard>::enqueue(const T &v) {
  bool done = false;

  size_t chead = _head.load(std::memory_order_acquire);
  while(!done) {
    size_t ctail = _tail.load(std::memory_order_acquire);
    elt celt;
    celt.binary = _buffer[chead % _size].load(std::memory_order_relaxed);
    size_t achead = _head.load(std::memory_order_acquire);
    if(achead != chead) {
      chead = achead;
      continue;
    }
    size_t nhead = chead + 1;
    if(nhead % _size == ctail % _size)
      return false;
    if(celt.split.val == guard) {
      // entry is empty
      elt nelt;
      nelt.split.id  = celt.split.id + 1;
      nelt.split.val = v;
      done = _buffer[chead % _size].compare_exchange_weak(celt.binary, nelt.binary);
      // done == true <=> sucessfully written entry
    }
    _head.compare_exchange_strong(chead, nhead);
    //    cas(&_head, chead, nhead, &chead);
  }

  return true;
}

template<typename T, unsigned int n, T guard>
void jflib::basic_circular_buffer<T,n,guard>::enqueue_no_check(const T &v) {
  bool done = false;

  size_t chead = _head.load(std::memory_order_acquire);
  while(!done) {
    elt celt;
    celt.binary = _buffer[chead % _size].load(std::memory_order_relaxed);
    size_t achead = _head.load(std::memory_order_acquire);
    if(achead != chead) {
      chead = achead;
      continue;
    }
    size_t nhead = chead + 1;
    if(celt.split.val == guard) {
      // entry is empty
      elt nelt;
      nelt.split.id  = celt.split.id + 1;
      nelt.split.val = v;
      //      done = cas(&_buffer[chead % _size], celt.binary, nelt.binary);
      done = _buffer[chead % _size].compare_exchange_weak(celt.binary, nelt.binary);
      // done == true <=> sucessfully written entry
    }
    // cas(&_head, chead, nhead, &chead);
    _head.compare_exchange_strong(chead, nhead);
  }
}

template<typename T, unsigned int n, T guard>
T jflib::basic_circular_buffer<T,n,guard>::dequeue() {
  bool done = false;
  elt res;

  size_t ctail = _tail.load(std::memory_order_acquire);
  while(!done) {
    bool dequeued = false;
    do {
      if(ctail % _size == _head.load(std::memory_order_acquire) % _size)
        return guard;
      size_t ntail = ctail + 1;
      //      dequeued = cas(&_tail, ctail, ntail, &ctail);
      dequeued = _tail.compare_exchange_weak(ctail, ntail);
    } while(!dequeued);

    res.binary = _buffer[ctail % _size].load(std::memory_order_relaxed);
    elt nres;
    nres.split.val = guard;
    while(true) {
      nres.split.id = res.split.id + 1;
      if(res.split.val == guard) {
        //        if(cas(&_buffer[ctail % _size], res.binary, nres.binary, &res.binary))
        if(_buffer[ctail % _size].compare_exchange_strong(res.binary, nres.binary))
          break;
      } else {
        //        done = cas(&_buffer[ctail % _size], res.binary, nres.binary);
        done = _buffer[ctail % _size].compare_exchange_strong(res.binary, nres.binary);
        break;
      }
    }
  }

  return res.split.val;
}
#endif /* __JELLYFISH_CIRCULAR_BUFFER_HPP__ */
