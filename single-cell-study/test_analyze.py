import unittest
import numpy as np
import pandas as pd
from scipy import sparse
from analyze import aggregate, fit


class PairedAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.obs=pd.DataFrame([{'donor':f'd{d}','condition':c,'cell_type':'B'}
            for d in range(4) for c in ['ctrl','stim'] for _ in range(3)])
        self.obs.index=[f'c{i}' for i in range(len(self.obs))]
        self.x=sparse.csr_matrix(np.arange(len(self.obs)*2).reshape(-1,2))

    def test_exact_sums_and_complete_pairs(self):
        counts,meta,audit=aggregate(self.x,self.obs,['g1','g2'],'B',3)
        np.testing.assert_array_equal(counts.iloc[0],np.asarray(self.x[:3].sum(0)).ravel())
        self.assertEqual(len(meta),8)
        keep=~((self.obs.donor=='d0')&(self.obs.condition=='stim'))
        counts,meta,audit=aggregate(self.x[keep.to_numpy()],self.obs[keep],['g1','g2'],'B',3)
        self.assertEqual(len(meta),6)
        self.assertFalse(audit.loc['d0','retained'])

    def test_reject_missing_and_fractional(self):
        bad=self.obs.copy();bad.loc['c0','donor']=None
        with self.assertRaises(ValueError):aggregate(self.x,bad,['g1','g2'],'B',3)
        with self.assertRaises(ValueError):aggregate(self.x*.5,self.obs,['g1','g2'],'B',3)

    def test_known_paired_effect(self):
        rng=np.random.default_rng(17)
        meta=pd.DataFrame([{'donor':f'd{d}','condition':c} for d in range(4) for c in ['ctrl','stim']])
        mu=np.full((8,200),100.0);mu[meta.condition=='stim',:20]=1600
        counts=pd.DataFrame(rng.negative_binomial(20,20/(20+mu)),index=meta.index,
                            columns=[f'g{i}' for i in range(200)])
        res,_=fit(counts,meta)
        self.assertGreater(res.loc[[f'g{i}' for i in range(20)],'log2FoldChange'].median(),2)
        self.assertTrue((res.padj<.05).any())


if __name__=='__main__':unittest.main()
