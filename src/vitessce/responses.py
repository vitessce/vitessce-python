import ujson
from starlette.responses import JSONResponse

# References:
# - https://www.starlette.io/responses/#custom-json-serialization
# - https://github.com/encode/starlette/releases/tag/0.14.1


class UJSONResponse(JSONResponse):
    def render(self, content):
        return ujson.dumps(content, ensure_ascii=False).encode("utf-8")
