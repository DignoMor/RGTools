
import sys

class Logger:
    def __init__(self, opath="stderr"):
        self._opath = opath
    
    def take_log(self, message):
        if self._opath == "stderr":
            sys.stderr.write(message)
            sys.stderr.write("\n")

        elif self._opath == "stdout":
            sys.stdout.write(message)
            sys.stdout.write("\n")

        else:
            with open(self._opath, "a") as f:
                f.write(message)
                f.write("\n")
