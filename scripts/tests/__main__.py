import os
import sys

import pytest

if __name__ == "__main__":
    
    test_dir = os.path.dirname(os.path.abspath(__file__)) # all test in the folder

    sys.exit(pytest.main([test_dir, "-v"]))