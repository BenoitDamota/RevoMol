import os
import sys

# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


from evomol.main import main

if __name__ == "__main__":
    main()
