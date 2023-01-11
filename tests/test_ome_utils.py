import unittest
from pathlib import Path
import zarr
import numpy as np

from vitessce.data_utils import (
    rgb_img_to_ome_zarr,
)


data_path = Path('tests/data')


class TestOmeUtils(unittest.TestCase):

    def setUp(self):
        self.img_arr = np.array([
            [
                [150, 149, 150, 150, 149, 150, 150, 150],
                [150, 150, 149, 149, 149, 149, 150, 150],
                [149, 149, 149, 150, 150, 149, 149, 150],
                [148, 148, 149, 149, 149, 149, 149, 149],
                [150, 149, 150, 149, 150, 149, 148, 149],
                [149, 149, 150, 149, 150, 150, 149, 149],
                [150, 149, 150, 150, 149, 149, 150, 148],
                [150, 149, 149, 150, 149, 149, 148, 147]
            ],
            [
                [153, 153, 153, 153, 152, 152, 152, 152],
                [153, 153, 153, 153, 153, 153, 152, 153],
                [153, 153, 152, 153, 153, 152, 152, 153],
                [153, 152, 152, 152, 153, 153, 153, 153],
                [153, 152, 152, 152, 153, 153, 153, 153],
                [152, 152, 153, 152, 153, 152, 153, 152],
                [153, 152, 152, 152, 153, 153, 153, 152],
                [152, 152, 152, 153, 153, 153, 153, 153]
            ],
            [
                [147, 145, 145, 145, 146, 146, 147, 148],
                [146, 147, 147, 148, 148, 147, 146, 147],
                [148, 147, 146, 147, 148, 148, 147, 148],
                [148, 147, 146, 147, 148, 148, 148, 148],
                [147, 147, 146, 146, 146, 145, 146, 147],
                [146, 148, 146, 146, 146, 146, 146, 148],
                [147, 146, 145, 146, 147, 146, 146, 147],
                [147, 147, 147, 147, 148, 147, 147, 147]
            ]
        ])

    def test_rgb_img_to_ome_zarr(self):
        img_arr = self.img_arr
        out_path = data_path / "rgb_out.ome.zarr"
        rgb_img_to_ome_zarr(img_arr, out_path, img_name="Test", axes="cyx", chunks=(1, 3, 3))

        z_root = zarr.open(out_path, mode="r")

        assert dict(z_root.attrs) == {
            'multiscales': [
                {
                    'axes': ['c', 'y', 'x'],
                    'datasets': [
                        {'path': '0'},
                        {'path': '1'},
                        {'path': '2'},
                        {'path': '3'},
                        {'path': '4'}
                    ],
                    'version': '0.3'
                }
            ],
            'omero': {
                'channels': [
                    {
                        'color': 'FF0000',
                        'label': 'R',
                        'window': {'end': 255, 'max': 255, 'min': 0, 'start': 0}
                    },
                    {
                        'color': '00FF00',
                        'label': 'G',
                        'window': {'end': 255, 'max': 255, 'min': 0, 'start': 0}
                    },
                    {
                        'color': '0000FF',
                        'label': 'B',
                        'window': {'end': 255, 'max': 255, 'min': 0, 'start': 0}
                    }
                ],
                'name': 'Test',
                'rdefs': {},
                'version': '0.3'
            }
        }
