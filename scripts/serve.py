"""
FactorForge local server — serves web/ static files + /api/optimize on a single port.

Usage:
    python server.py          # http://localhost:8080
    PORT=9000 python server.py
"""

from __future__ import annotations

import os
import sys
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

ROOT = Path(__file__).parent
sys.path.insert(0, str(ROOT / "src"))
os.chdir(ROOT / "web")

sys.path.insert(0, str(ROOT))
from api.optimize import handler as OptimizeAPIHandler  # noqa: E402


class FactorForgeHandler(OptimizeAPIHandler, SimpleHTTPRequestHandler):
    """Serves web/ for GET and routes /api/* to the optimize handler.

    Inherits from OptimizeAPIHandler (not just delegating to it via an
    unbound method call) so that self.validate_host / self.send_error_response
    and any other instance methods the optimize handler relies on resolve
    correctly. A previous version called `handler.do_POST(self)` with a
    foreign `self`, which broke as soon as optimize.py's do_POST started
    calling additional self.* methods not present on this class.
    """

    def do_POST(self) -> None:
        if self.path == "/api/optimize":
            OptimizeAPIHandler.do_POST(self)
        else:
            self.send_error(404, "Not found")

    def do_GET(self) -> None:
        if self.path.startswith("/api/"):
            OptimizeAPIHandler.do_GET(self)
        else:
            SimpleHTTPRequestHandler.do_GET(self)

    def log_message(self, fmt: str, *args: object) -> None:
        print(f"[FactorForge] {fmt % args}")


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8080))
    server = HTTPServer(("0.0.0.0", port), FactorForgeHandler)
    print(f"FactorForge running at http://localhost:{port}")
    print("Press Ctrl+C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
