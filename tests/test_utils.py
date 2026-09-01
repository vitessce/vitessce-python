import unittest
from urllib.parse import unquote

from vitessce.utils import make_ids_csv_data_url, make_colors_csv_data_url


class TestCsvDataUrlUtils(unittest.TestCase):

    def test_make_ids_csv_data_url(self):
        url = make_ids_csv_data_url([612, 3351, 4328])
        self.assertTrue(url.startswith("data:text/csv,"))
        decoded = unquote(url[len("data:text/csv,"):])
        self.assertEqual(decoded, "id\r\n612\r\n3351\r\n4328\r\n")

    def test_make_ids_csv_data_url_empty(self):
        url = make_ids_csv_data_url([])
        decoded = unquote(url[len("data:text/csv,"):])
        self.assertEqual(decoded, "id\r\n")

    def test_make_colors_csv_data_url(self):
        url = make_colors_csv_data_url({
            612: "#d74242",
            3351: "#b9d742",
        })
        self.assertTrue(url.startswith("data:text/csv,"))
        decoded = unquote(url[len("data:text/csv,"):])
        self.assertEqual(decoded, "id,color\r\n612,#d74242\r\n3351,#b9d742\r\n")


if __name__ == "__main__":
    unittest.main()
