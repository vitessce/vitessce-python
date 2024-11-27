
from starlette.routing import Route
from starlette.responses import StreamingResponse
from pathlib import Path

# Adapted from https://gist.github.com/tombulled/712fd8e19ed0618c5f9f7d5f5f543782


def ranged(file, start=0, end=None, block_size=65535):
    consumed = 0
    file.seek(start)

    while True:
        data_length = (
            min(block_size, end - start - consumed)
            if end else block_size
        )
        if data_length <= 0:
            break
        data = file.read(data_length)
        if not data:
            break
        consumed += data_length
        yield data

    if hasattr(file, 'close'):
        file.close()


def range_repsonse(request, file_path):
    path = Path(file_path)
    file = path.open('rb')
    file_size = path.stat().st_size
    content_range = request.headers.get('range')
    content_length = file_size
    status_code = 200
    headers = {}

    if content_range is not None:
        content_range = content_range.strip().lower()
        content_ranges = content_range.split('=')[-1]
        range_start, range_end, * \
            _ = map(str.strip, (content_ranges + '-').split('-'))
        range_start = max(0, int(range_start)) if range_start else 0
        range_end = min(file_size - 1, int(range_end)
                        ) if range_end else file_size - 1
        content_length = (range_end - range_start) + 1
        file = ranged(file, start=range_start, end=range_end + 1)
        status_code = 206
        headers['Content-Range'] = f'bytes {range_start}-{range_end}/{file_size}'

    response = StreamingResponse(
        file,
        media_type='tiff',
        status_code=status_code,
    )
    response.headers.update({
        'Accept-Ranges': 'bytes',
        'Content-Length': str(content_length),
        **headers,
    })
    return response


class JsonRoute(Route):
    def __init__(self, path, endpoint, data_json):
        super().__init__(path, endpoint)
        self.data_json = data_json


class FileRoute(Route):
    def __init__(self, path, endpoint, file_path):
        super().__init__(path, endpoint)
        self.file_path = file_path
