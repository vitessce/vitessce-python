import json
import os
from os.path import join
import tempfile

from starlette.routing import Route, Mount
from starlette.staticfiles import StaticFiles

from .routes import create_obj_routes
from .wrappers import JsonRoute


def upload_to_s3(config, s3, bucket_name, prefix=''):
    """
    :param config: The Vitessce view config to upload.
    :type config: VitessceConfig
    :param s3: A boto3 S3 resource object with permission to upload to the specified bucket.
    :type s3: boto3.resource
    :param str bucket_name: The name of the bucket to which to upload.
    :param str prefix: The prefix path for the bucket keys (think subdirectory).
    """
    
    base_url = f"https://{bucket_name}.s3.amazonaws.com" + ("/" + prefix if len(prefix) > 0 else "")
    bucket = s3.Bucket(bucket_name)
    
    routes = []
    def on_obj(obj, dataset_uid, obj_i):
        obj_file_defs, obj_routes = create_obj_routes(obj, base_url, dataset_uid, obj_i)
        for obj_route in obj_routes:
            routes.append(obj_route)
        return obj_file_defs
    config_dict = config.to_dict(on_obj=on_obj)

    for route in routes:
        route_path = route.path[1:]
        key = (prefix + "/" if len(prefix) > 0 else "") + route_path

        print(f"Uploading {bucket_name}:{key}")
        
        if type(route) == JsonRoute:
            data_json = route.data_json
            with tempfile.TemporaryFile() as f:
                f.write(json.dumps(data_json).encode())
                f.write('\n'.encode())
                f.flush()
                f.seek(0, 0)
                bucket.put_object(Key=key, Body=f)
        elif type(route) == Mount:
            route_app = route.app
            if type(route_app) == StaticFiles:
                static_dir = route_app.directory

                for root, dirs, files in os.walk(static_dir):
                    for filename in files:
                        file_key = key + join(root, filename)[len(static_dir):]
                        bucket.upload_file(join(root, filename), file_key)
    
    return config_dict