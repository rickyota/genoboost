""" run by `python -m unittest discover -s ./projects/genetics_py`
 not sure how to run all projects at once
"""

import unittest
from io import StringIO
import pandas as pd


from src import plink


class TestDataset(unittest.TestCase):
    def test_plink(self):
        fam = 'a 1 0 0 1 2\nb 2 0 0 2 1\n'
        #fam = 'a 1 0 0 1 1 1\nb 2 0 0 2 2\n'
        fam = plink.load_fam(StringIO(fam))

        cols = 'fid,iid,sex,phe'.split(',')
        df_ans = pd.DataFrame([['a', '1', 1, 2], ['b', '2', 2, 1]], columns=cols)
        pd.testing.assert_frame_equal(fam[cols], df_ans)


if __name__ == '__main__':
    unittest.main()
