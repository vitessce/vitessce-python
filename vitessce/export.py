import json
import os
from os.path import join
from shutil import copyfile

from starlette.routing import Mount
from starlette.staticfiles import StaticFiles

from .routes import JsonRoute, FileRoute


def export_to_s3(config, s3, bucket_name, prefix=''):
    """
    :param config: The Vitessce view config to export to S3.
    :type config: VitessceConfig
    :param s3: A boto3 S3 resource object with permission to upload to the specified bucket.
    :type s3: boto3.resource
    :param str bucket_name: The name of the bucket to which to upload.
    :param str prefix: The prefix path for the bucket keys (think subdirectory).

    :returns: The config as a dict, with S3 urls filled in.
    :rtype: dict
    """

    base_url = f"https://{bucket_name}.s3.amazonaws.com" + \
        ("/" + prefix if len(prefix) > 0 else "")
    bucket = s3.Bucket(bucket_name)
    config_dict = config.to_dict(base_url=base_url)
    routes = config.get_routes()
    uploaded_routes = []
    for route in routes:
        route_path = route.path[1:]
        key = (prefix + "/" if len(prefix) > 0 else "") + route_path

        print(f"Uploading {bucket_name}:{key}")

        if isinstance(route, JsonRoute):
            if route not in uploaded_routes:
                data_json = route.data_json
                bucket.put_object(Key=key, Body=json.dumps(data_json).encode())
                uploaded_routes.append(route)
        elif isinstance(route, FileRoute):
            if route not in uploaded_routes:
                local_file_path = route.file_path
                s3.meta.client.upload_file(local_file_path, bucket_name, key)
                uploaded_routes.append(route)
        elif isinstance(route, Mount):
            route_app = route.app
            if isinstance(route_app, StaticFiles):
                if route not in uploaded_routes:
                    uploaded_routes.append(route)
                    static_dir = route_app.directory
                    for root, dirs, files in os.walk(static_dir):
                        for filename in files:
                            file_key = key + \
                                join(root, filename)[len(static_dir):]
                            bucket.upload_file(join(root, filename), file_key)

    return config_dict


def export_to_files(config, base_url, out_dir='.'):
    """
    :param config: The Vitessce view config to export to files.
    :type config: VitessceConfig
    :param str out_dir: The path to the output directory. By default, the current directory.
    :param str base_url: The URL on which the files will be served.

    :returns: The config as a dict, with urls filled in.
    :rtype: dict
    """

    config_dict = config.to_dict(base_url=base_url)
    routes = config.get_routes()
    for route in routes:
        route_path = route.path[1:]
        out_path = join(out_dir, route_path)

        if isinstance(route, JsonRoute):
            data_json = route.data_json
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            with open(out_path, 'w') as f:
                json.dump(data_json, f)
        elif isinstance(route, FileRoute):
            local_file_path = route.file_path
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            copyfile(local_file_path, out_path)
        elif isinstance(route, Mount):
            route_app = route.app
            if isinstance(route_app, StaticFiles):
                static_dir = route_app.directory

                for root, dirs, files in os.walk(static_dir):
                    for filename in files:
                        file_key = out_path + \
                            join(root, filename)[len(static_dir):]
                        os.makedirs(os.path.dirname(file_key), exist_ok=True)
                        copyfile(join(root, filename), file_key)

    return config_dict
