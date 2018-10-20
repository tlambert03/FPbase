import json
from .skbio_protein import SkbSequence
from .util import protein_weight, slugify
from .align import nw_align
from .mutations import Mutations, detect_offset

try:
    import requests
except ImportError:
    print('Could not import requests. Cannot pull sequences from fpbase')


class FPSeq(SkbSequence):

    def __init__(self, sequence, position_lables=None, **kwargs):
        if isinstance(position_lables, dict):
            pos_labels = list()
            i = 1
            for n in range(i, len(sequence) + i):
                if n in position_lables:
                    pos_labels.append(str(position_lables[n]))
                else:
                    pos_labels.append(str(i))
                    i += 1
            kwargs['positional_metadata'] = {'position_lables': pos_labels}
        super().__init__(sequence, **kwargs)

    @property
    def weight(self):
        try:
            return protein_weight(str(self)) / 1000
        except ValueError:
            pass

    def align_to(self, other, **kwargs):
        return nw_align(str(self), str(other), **kwargs)

    def mutations_to(self, other, **kwargs):
        return self.align_to(other, **kwargs).mutations

    def mutate(self, mutations, zeroindex=1, err_on_shift=False):
        if isinstance(mutations, str):
            mutations = Mutations(mutations)
        offset = detect_offset(self, mutations)
        if offset != 0:
            if offset is None or err_on_shift:
                raise ValueError('Requested mutation numbers inconsistent '
                                 'with current sequence frame')
            else:
                print('WARNING: mutations shifted {} position{} relative to sequence frame'
                      .format(offset, 's' if abs(offset) > 1 else ''))
        out = list(str(self))
        for before, index, after in mutations:
                out[index - zeroindex + offset] = after
        return FPSeq("".join(out))

    @classmethod
    def from_fpbase(cls, slug):
        url = 'https://www.fpbase.org/api/{}/?format=json'.format(slugify(slug))
        response = requests.get(url)
        return FPSeq(json.loads(response.content).get('seq'))
