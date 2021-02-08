from pathlib import Path

import rebound

Path("librebound.so").unlink(missing_ok=True)
Path("rebound.h").unlink(missing_ok=True)

module_dir = Path(rebound.__file__).parent
rebound_h = module_dir / "rebound.h"
librebound = Path(rebound.__libpath__)

print(librebound)
print(rebound_h)

Path("librebound.so").symlink_to(librebound)
Path("rebound.h").symlink_to(rebound_h)
