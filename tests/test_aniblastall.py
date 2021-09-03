from pathlib import Path

from pyani import aniblastall


# Test get_version()
# Test case 1: there is no executable
def test_get_version_missing_exe(executable_missing):
    """Test behaviour when there is no file at the specified executable location."""
    test_file_1 = Path("/non/existent/blastall")
    assert (
        aniblastall.get_version(test_file_1)
        == f"No blastall executable at {test_file_1}"
    )


# Test case 2: there is a file, but it is not executable
def test_get_version_not_executable(executable_not_executable):
    """Test behaviour when the file at the executable location is not executable."""
    test_file_2 = Path("/non/executable/blastall")
    assert (
        aniblastall.get_version(test_file_2)
        == f"blastall exists at {test_file_2} but not executable"
    )


# Test case 3: there is an executable file, but the version can't be retrieved
def test_get_version_no_version(executable_without_version):
    """Test behaviour when the version for the executable can not be retrieved."""
    test_file_3 = Path("/missing/version/blastall")
    assert (
        aniblastall.get_version(test_file_3)
        == f"blastall exists at {test_file_3} but could not retrieve version"
    )


# Test case 4: there is an executable file, but it will not run on the OS
def test_get_version_os_incompatible(executable_incompatible_with_os):
    """Test behaviour when the program can't run on the operating system.
    This will happen with newer versions of MacOS."""
    test_file_4 = Path("/os/incompatible/blastall")
    assert (
        aniblastall.get_version(test_file_4)
        == f"blastall exists at {test_file_4} but could not be executed"
    )
