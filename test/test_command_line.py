import unittest
from unittest.mock import patch

from asari import __version__
from asari.command_line import main


class TestCommandLine(unittest.TestCase):

    def test_no_subcommand_exits(self):
        # subcommand is required; argparse exits with an error
        with patch('sys.argv', ['asari']), patch('sys.stderr'):
            with self.assertRaises(SystemExit) as cm:
                main()
        self.assertNotEqual(cm.exception.code, 0)

    def test_version(self):
        with patch('sys.argv', ['asari', '-v']), patch('sys.stdout') as out:
            with self.assertRaises(SystemExit) as cm:
                main()
        self.assertEqual(cm.exception.code, 0)
        self.assertIn(__version__, ''.join(c.args[0] for c in out.write.call_args_list))


if __name__ == '__main__':
    unittest.main()
