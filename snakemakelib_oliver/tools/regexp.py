import os
from snakemakelib.bio.ngs.regexp import RegexpDict


class SampleRegexp(RegexpDict):
    _required_keys = ['SM']
    _group_keys = ['PU']
    _extra_keys = ['PATH', 'FILE']

    def __init__(self, regexp=None, *args, **kwargs):
        super(SampleRegexp, self).__init__(regexp, *args, **kwargs)

    def _post_process_keys(self, m):
        self['PATH'] = os.path.dirname(m.string)


class RunRegexp(RegexpDict):
    _group_keys = ['SM', 'PU', 'DT']
    _extra_keys = ['PATH', 'FILE']

    def __init__(self, regexp=None, *args, **kwargs):
        super(RunRegexp, self).__init__(regexp, *args, **kwargs)

    def _post_process_keys(self, m):
        self['PATH'] = os.path.dirname(m.string)

