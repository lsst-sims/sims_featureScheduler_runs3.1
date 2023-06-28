import numpy as np
import rubin_sim.maf as maf
import os
import argparse
import rubin_sim.maf.db as db


class GapsMetric(maf.BaseMetric):
    """Compute the number of the time gaps between observations that occur in a given time range.

    Parameters
    ----------
    times_col : `str`, opt
        The column name for the exposure times.  Values assumed to be in days.
        Default observationStartMJD.
    time_scale : `float` (2/24)
        Time scale to see how well it is sampled (days).
   

    Returns
    -------
    """

    def __init__(
        self,
        times_col="observationStartMJD",
        time_scale=2. / 24,
        units="N",
        **kwargs,
    ):
        self.times_col = times_col
        # divide by two so we bin at the Nyquist frequency
        self.bin_size = time_scale / 2.
        super().__init__(
            col=[self.times_col], metric_dtype="float", units=units, **kwargs
        )

    def run(self, data_slice, slice_point=None):
        if data_slice.size < 2:
            return self.badval
        times = np.sort(data_slice[self.times_col])
        bins = np.arange(times.min()-self.bin_size, times.max()+self.bin_size, self.bin_size)
        vals, _be = np.histogram(times, bins)
        mult = vals[2:] * vals[0:-2]
        result = np.size(np.where(mult > 0)[0])
        
        return result


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=str, default=None)
    args = parser.parse_args()

    filename = args.db
    run_name = os.path.basename(filename).replace('.db', '')
    out_dir = run_name + "_sci"
    results_db = db.ResultsDb(out_dir=run_name + "_sci")

    scales = [2./24, 5./25, 1., 2, 3, 4]

    bundle_list = []
    for time_scale in scales:
        for filtername in 'ugrizy':
            metric = GapsMetric(time_scale=time_scale, metric_name='gaps %s %i' % (filtername, time_scale*24))
            slicer = maf.HealpixSlicer()
            summary_sats = [maf.MedianMetric(), maf.SumMetric()]
            plot_dict = {'percentile_clip': 95.}
            sql = "filter='%s' and note not like '%%neo%%'" % filtername

            bundle_list.append(maf.MetricBundle(metric, slicer, sql, plot_dict=plot_dict,
                                                summary_metrics=summary_sats, run_name=run_name))

    mbg = maf.MetricBundleGroup(bundle_list, filename, out_dir=out_dir, results_db=results_db)
    mbg.run_all()
    mbg.plot_all()
