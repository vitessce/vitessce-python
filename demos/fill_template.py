import os
from os.path import join
from jinja2 import Environment, FileSystemLoader
import argparse


def render_json(dir_name, version_str, url_type, port):
    template_dir = join(os.path.dirname(os.path.abspath(__file__)), dir_name)
    env = Environment(loader=FileSystemLoader(template_dir))

    template = env.get_template('vitessce.template.json')

    BASE_URL = {
        'local': f'http://localhost:{port}/{dir_name}/data/processed',
        'remote': f'https://s3.amazonaws.com/vitessce-data/{version_str}/main/{dir_name}'
    }
    BASE_URL_GCP = {
        'local': f'http://localhost:{port}/{dir_name}/data/processed',
        'remote': f'https://vitessce-data.storage.googleapis.com/{version_str}/main/{dir_name}'
    }

    out = template.render(
        base_url=BASE_URL[url_type],
        base_url_gcp=BASE_URL_GCP[url_type]
    )
    print(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d',
        '--dir_name',
        type=str,
        required=True,
        help='The demo subdirectory'
    )
    parser.add_argument(
        '-t',
        '--url_type',
        type=str,
        required=True,
        help='The URL type, either local or remote'
    )
    parser.add_argument(
        '-v',
        '--version',
        type=str,
        required=False,
        help='The version'
    )
    parser.add_argument(
        '-p',
        '--port',
        type=str,
        default='8000',
        required=False,
        help='The port to use for local URLs'
    )
    args = parser.parse_args()
    if args.url_type == 'remote' and args.version is None:
        raise ValueError('The --version argument must be provided when --url_type == "remote"')
    render_json(
        args.dir_name,
        args.version,
        args.url_type,
        args.port
    )
