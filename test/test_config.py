import json
import unittest

from vitessce import VitessceConfig

class TestConfig(unittest.TestCase):

    def test_config_creation(self):
        vc = VitessceConfig()
        vc_dict = vc.to_dict()
        vc_json = json.dumps(vc_dict)

        self.assertEqual(vc_dict, {
            "version": "1.0.0",
            "name": "",
            "description": "",
            "datasets": [],
            "coordinationSpace": {
                "dataset": {},
                "embeddingType": {},
            },
            "layout": [],
            "initStrategy": "auto"
        })
    
    def test_config_add_dataset(self):
        vc = VitessceConfig()
        my_dataset = vc.add_dataset(name='My Dataset')

        vc_dict = vc.to_dict()
        vc_json = json.dumps(vc_dict)

        self.assertEqual(len(vc_dict["datasets"]), 1)
        self.assertEqual(len(vc_dict["coordinationSpace"]["dataset"]), 1)

        self.assertEqual(vc_dict, {
            "version": "1.0.0",
            "name": "",
            "description": "",
            "datasets": [
                { 'uid': 'A', 'name': 'My Dataset', 'files': [] }
            ],
            'coordinationSpace': {
                'dataset': {'A': 'A'},
                'embeddingType': {}
            },
            "layout": [],
            "initStrategy": "auto"
        })

    def test_load_config(self):

        vc = VitessceConfig(config={
            "version": "1.0.0",
            "name": "Test name",
            "description": "Test description",
            "datasets": [
                { 'uid': 'A', 'name': 'My First Dataset', 'files': [
                    {
                        'url': 'http://cells.json',
                        'type': 'cells',
                        'fileType': 'cells.json'
                    }
                ] }
            ],
            'coordinationSpace': {
                'dataset': {'A': 'A'},
                'embeddingType': {}
            },
            "layout": [],
            "initStrategy": "auto"
        })

        my_second_dataset = vc.add_dataset(name='My Second Dataset')
        vc_dict = vc.to_dict()
        vc_json = json.dumps(vc_dict)

        self.assertEqual(len(vc_dict["datasets"]), 2)
        self.assertEqual(vc_dict, {
            "version": "1.0.0",
            "name": "Test name",
            "description": "Test description",
            "datasets": [
                { 'uid': 'A', 'name': 'My First Dataset', 'files': [{
                        'url': 'http://cells.json',
                        'type': 'cells',
                        'fileType': 'cells.json'
                    }] },
                { 'uid': 'B', 'name': 'My Second Dataset', 'files': [] }
            ],
            'coordinationSpace': {
                'dataset': {
                    'A': 'A',
                    'B': 'B'
                },
                'embeddingType': {}
            },
            "layout": [],
            "initStrategy": "auto"
        })

