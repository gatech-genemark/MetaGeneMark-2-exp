import logging
import pandas as pd
from typing import *

from sbsp_io.assembly_summary import read_assembly_summary_into_dataframe

log = logging.getLogger(__name__)

class AssemblySummary(pd.DataFrame):

    @property
    def _constructor(self):
        """Used so that this type is always returned by pandas actions"""
        return AssemblySummary

    @classmethod
    def init_from_file(cls, pf_assembly_summary):
        # type: (str) -> AssemblySummary
        df = read_assembly_summary_into_dataframe(pf_assembly_summary)
        return cls(df)

    def write(self, pf_output):
        # type: (str) -> None

        df = self.rename(columns={"assembly_accession": "# assembly_accession"})
        df.to_csv(pf_output, sep="\t", index=False)

