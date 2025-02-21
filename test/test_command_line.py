import unittest
from unittest.mock import patch
from asari.command_line import main

class TestCommandLine(unittest.TestCase):

    @patch('asari.command_line.main', autospec=True)
    def test_main_called(self, mock_main):
        # Call the main function and expect SystemExit
        with self.assertRaises(SystemExit):
            main([])

if __name__ == '__main__':
    unittest.main()