import unittest
from unittest import mock

from fragutils.utils.network_utils import get_driver


@mock.patch("neo4j.GraphDatabase.driver")
class GetDriverTest(unittest.TestCase):
    def test_bare_hostname_uses_plain_bolt(self, mock_driver):
        get_driver("graph.graph-a.svc", "neo4j/secret")

        mock_driver.assert_called_once_with(
            "bolt://graph.graph-a.svc:7687", auth=("neo4j", "secret")
        )

    def test_bare_hostname_without_auth(self, mock_driver):
        get_driver("graph.graph-a.svc", "none")

        mock_driver.assert_called_once_with("bolt://graph.graph-a.svc:7687")

    def test_uri_with_secure_scheme_is_used_unchanged(self, mock_driver):
        uri = "bolt+s://graph-x.xchem-dev.diamond.ac.uk:7687"

        get_driver(uri, "neo4j/secret")

        mock_driver.assert_called_once_with(uri, auth=("neo4j", "secret"))

    def test_uri_without_auth_is_used_unchanged(self, mock_driver):
        uri = "neo4j+s://graph-x.xchem-dev.diamond.ac.uk:7687"

        get_driver(uri, "none")

        mock_driver.assert_called_once_with(uri)


if __name__ == "__main__":
    unittest.main()
