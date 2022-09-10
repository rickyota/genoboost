""" run by `python -m unittest discover -s ./projects/genetics_py`
 not sure how to run all projects at once
"""

import unittest
import pandas as pd


from src import dataset


class TestDataset(unittest.TestCase):
    def test_split_obs_ts(self):
        n = 100
        test_prop = 0.2

        sample_id = [str(i) for i in range(n)]
        sample_id = {'fid': sample_id, 'iid': sample_id}
        sample_id = pd.DataFrame(sample_id)

        obs_id, test_id = dataset.split_obs_ts(sample_id, test_prop)
        self.assertEqual(len(test_id), int(n * test_prop))

        intersection = pd.merge(obs_id, test_id, how='inner')
        self.assertEqual(len(intersection), 0)

        union = pd.concat([obs_id, test_id], join='outer')
        union = union.sort_index()
        pd.testing.assert_frame_equal(union, sample_id)

    def test_split_tr_va_cv(self):
        n = 100
        cv_n = 5
        sample_id = [str(i) for i in range(n)]
        sample_id = {'fid': sample_id, 'iid': sample_id}
        sample_id = pd.DataFrame(sample_id)

        ids = dataset.split_tr_va_cv(sample_id, cv_n=cv_n)

        for (_, (tr, va)) in enumerate(ids):
            self.assertEqual(len(va), int(n / cv_n))

            intersection = pd.merge(tr, va, how='inner')
            self.assertEqual(len(intersection), 0)

            union = pd.concat([tr, va], join='outer')
            union = union.sort_index()
            pd.testing.assert_frame_equal(union, sample_id)


if __name__ == '__main__':
    unittest.main()
