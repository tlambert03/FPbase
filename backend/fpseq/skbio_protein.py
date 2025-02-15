# The vast majority of this code is borrowed form the SciKit bio package
# with the goal of minimizing dependencies... I'm copying just what I need

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice, this
#   list of conditions and the following disclaimer in the documentation and/or
#   other materials provided with the distribution.
#
# * Neither the names scikit-bio, skbio, or biocore nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ----------------------------------------------------------------------------


import numbers

import numpy as np

from .util import chunk_string


class classproperty(property):
    """Decorator for class-level properties.
    borrowed from scikit-bio
    https://github.com/biocore/scikit-bio/blob/master/skbio/util/_decorator.py#L304
    """

    def __init__(self, func):
        name = func.__name__
        doc = func.__doc__
        super().__init__(classmethod(func))
        self.__name__ = name
        self.__doc__ = doc

    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()

    def __set__(self, obj, value):
        raise AttributeError("can't set attribute")


class SkbSequence:
    """trying not to import the full scikit-bio package... this is a minimal
    rebuild of the Scikit-Bio Grammared Sequence Class.  Please see sci-kit
    bio package for full Sequence class and documentation!!
    https://github.com/biocore/scikit-bio
    """

    _number_of_extended_ascii_codes = 256
    __validation_mask = None
    __degenerate_codes = None
    __definite_char_codes = None
    __gap_codes = None

    _bytes: np.ndarray

    def __init__(self, sequence, validate=True):
        if isinstance(sequence, np.ndarray):
            if sequence.dtype == np.uint8:
                self._set_bytes_contiguous(sequence)
            elif sequence.dtype == "|S1":
                sequence = sequence.view(np.uint8)
                if sequence.shape == ():
                    sequence = np.array([sequence], dtype=np.uint8)
                self._set_bytes_contiguous(sequence)
            else:
                raise TypeError(
                    "Can only create sequence from numpy.ndarray of dtype "
                    f"np.uint8 or '|S1'. Invalid dtype: {sequence.dtype}"
                )
        elif isinstance(sequence, SkbSequence):
            sequence._assert_can_cast_to(type(self))
            sequence = sequence._bytes
            self._owns_bytes = False
            self._set_bytes(sequence)
        else:
            if isinstance(sequence, str):
                sequence = sequence.replace(" ", "").replace("\n", "")
                sequence = sequence.encode("ascii")

            s = np.frombuffer(sequence, dtype=np.uint8)
            if isinstance(sequence, np.generic) and len(s) != 1:
                raise TypeError(f"Can cannot create a sequence with {type(sequence).__name__!r}")

            sequence = s
            self._owns_bytes = True
            self._set_bytes(sequence)

        if validate:
            self._validate()

    def _set_bytes_contiguous(self, sequence):
        """Munge the sequence data into a numpy array of dtype uint8."""
        if not sequence.flags["C_CONTIGUOUS"]:
            # https://github.com/numpy/numpy/issues/5716
            sequence = np.ascontiguousarray(sequence)
            self._owns_bytes = True
        else:
            self._owns_bytes = False
        self._set_bytes(sequence)

    def _set_bytes(self, sequence: np.ndarray):
        sequence.flags.writeable = False
        self._bytes = sequence

    @classmethod
    def _assert_can_cast_to(cls, target):
        if not (issubclass(cls, target) or issubclass(target, cls)):
            raise TypeError(f"Cannot cast {cls.__name__!r} as {target.__name__!r}.")

    def __repr__(self):
        return "Protein\n" + "-" * 54 + "\n" + "\n".join(chunk_string(str(self), 10, 55))

    def __str__(self):
        return self._bytes.tobytes().decode("ascii")

    def _munge_to_sequence(self, other, method):
        if isinstance(other, SkbSequence):
            if type(other) is not type(self):
                raise TypeError(
                    f"Cannot use {self.__class__.__name__} and {other.__class__.__name__} together with `{method}`"
                )
            else:
                return other

        return SkbSequence(other)

    def _munge_to_bytestring(self, other, method):
        if isinstance(other, bytes):
            return other
        elif isinstance(other, str):
            return other.encode("ascii")
        else:
            return str(self._munge_to_sequence(other, method))

    def __contains__(self, subsequence):
        """Determine if a subsequence is contained in this sequence."""
        return self._munge_to_bytestring(subsequence, "in") in str(self)

    def __eq__(self, other):
        """Determine if this sequence is equal to another."""
        return bool(str(self) == str(other))

    def __ne__(self, other):
        """Determine if this sequence is not equal to another."""
        return not (self == other)

    def __len__(self):
        """Return the number of characters in this sequence."""
        return self._bytes.size

    def __bool__(self):
        """Returns truth value (truthiness) of sequence."""
        return len(self) > 0

    def __iter__(self):
        """Iterate over positions in this sequence."""
        for i in range(len(self)):
            yield self[i]

    def __getitem__(self, indexable):
        """Slice this sequence."""
        if not isinstance(indexable, np.ndarray) and (
            (not isinstance(indexable, str)) and hasattr(indexable, "__iter__")
        ):
            indexable_ = indexable
            indexable = np.asarray(indexable)

            if indexable.dtype == object:
                indexable = list(indexable_)  # TODO: Don't blow out memory

                if not indexable:
                    # indexing with an empty list, so convert to ndarray and
                    # fall through to ndarray slicing below
                    indexable = np.asarray(indexable)
                else:
                    seq = np.concatenate(list(_slices_from_iter(self._bytes, indexable)))
                    # index = _as_slice_if_single_index(indexable)

                    # positional_metadata = None
                    # if self.has_positional_metadata():
                    #     pos_md_slices = list(_slices_from_iter(
                    #                          self.positional_metadata, index))
                    #     positional_metadata = pd.concat(pos_md_slices)

                    # metadata = None
                    # if self.has_metadata():
                    #    metadata = self.metadata

                    return self._constructor(sequence=seq)
                    # metadata=metadata,
                    # positional_metadata=positional_metadata)

        elif isinstance(indexable, str | bool):
            raise IndexError(f"Cannot index with {type(indexable).__name__} type: {indexable!r}")

        if isinstance(indexable, np.ndarray) and indexable.dtype == bool and len(indexable) != len(self):
            raise IndexError(
                f"An boolean vector index must be the same length as the sequence ({len(self)}, not {len(indexable)})."
            )

        if isinstance(indexable, np.ndarray) and indexable.size == 0:
            # convert an empty ndarray to a supported dtype for slicing a numpy
            # array
            indexable = indexable.astype(int)

        seq = self._bytes[indexable]
        # positional_metadata = self._slice_positional_metadata(indexable)

        # metadata = None
        # if self.has_metadata():
        #     metadata = self.metadata

        return self.__class__(sequence=seq)

    @classproperty
    def _validation_mask(cls):
        # TODO These masks could be defined (as literals) on each concrete
        # object. For now, memoize!
        if cls.__validation_mask is None:
            cls.__validation_mask = np.invert(
                np.bincount(
                    np.frombuffer(("".join(cls.alphabet)).encode(), dtype=np.uint8),
                    minlength=cls._number_of_extended_ascii_codes,
                ).astype(bool)
            )
        return cls.__validation_mask

    @classproperty
    def _degenerate_codes(cls):
        if cls.__degenerate_codes is None:
            degens = cls.degenerate_chars
            cls.__degenerate_codes = np.asarray([ord(d) for d in degens])
        return cls.__degenerate_codes

    @classproperty
    def _definite_char_codes(cls):
        if cls.__definite_char_codes is None:
            definite_chars = cls.definite_chars
            cls.__definite_char_codes = np.asarray([ord(d) for d in definite_chars])
        return cls.__definite_char_codes

    @classproperty
    def _gap_codes(cls):
        if cls.__gap_codes is None:
            gaps = cls.gap_chars
            cls.__gap_codes = np.asarray([ord(g) for g in gaps])
        return cls.__gap_codes

    def _validate(self):
        """https://github.com/biocore/scikit-bio/blob/0.5.4/skbio/sequence/_grammared_sequence.py#L340"""
        invalid_characters = (
            np.bincount(self._bytes, minlength=self._number_of_extended_ascii_codes) * self._validation_mask
        )
        if np.any(invalid_characters):
            bad = list(np.where(invalid_characters > 0)[0].astype(np.uint8).view("|S1"))
            raise ValueError(
                "Invalid character{} in sequence: {!r}. \n"
                "Valid characters: {!r}\n"
                "Note: Use `lowercase` if your sequence contains lowercase "
                "characters not in the sequence's alphabet.".format(
                    "s" if len(bad) > 1 else "",
                    [str(b.tostring().decode("ascii")) for b in bad] if len(bad) > 1 else bad[0],
                    list(self.alphabet),
                )
            )

    @classproperty
    def alphabet(cls):
        return cls.degenerate_chars | cls.definite_chars | cls.gap_chars | cls.stop_chars

    @classproperty
    def degenerate_map(cls):
        return {"B": set("DN"), "Z": set("EQ"), "X": set("ACDEFGHIKLMNPQRSTVWY")}

    @classproperty
    def degenerate_chars(cls):
        return set(cls.degenerate_map)

    @classproperty
    def definite_chars(cls):
        return set("ACDEFGHIKLMNPQRSTVWY")

    @classproperty
    def gap_chars(cls):
        return set("-.")

    @classproperty
    def _stop_codes(cls):
        if cls.__stop_codes is None:
            stops = cls.stop_chars
            cls.__stop_codes = np.asarray([ord(s) for s in stops])
        return cls.__stop_codes

    @classproperty
    def stop_chars(cls):
        return set("*")

    def gaps(self):
        return np.in1d(self._bytes, self._gap_codes)

    def has_gaps(self):
        return bool(self.gaps().any())

    def degenerates(self):
        return np.in1d(self._bytes, self._degenerate_codes)

    def has_degenerates(self):
        return bool(self.degenerates().any())

    def stops(self):
        return np.in1d(self._bytes, self._stop_codes)

    def has_stops(self):
        return bool(self.stops().any())

    def degap(self):
        return self[np.invert(self.gaps())]


def _single_index_to_slice(start_index):
    end_index = None if start_index == -1 else start_index + 1
    return slice(start_index, end_index)


def _is_single_index(index):
    return isinstance(index, numbers.Integral) and not isinstance(index, bool)


def _as_slice_if_single_index(indexable):
    if _is_single_index(indexable):
        return _single_index_to_slice(indexable)
    else:
        return indexable


def _slices_from_iter(array, indexables):
    for i in indexables:
        if isinstance(i, slice):
            pass
        elif _is_single_index(i):
            i = _single_index_to_slice(i)
        else:
            raise IndexError(f"Cannot slice sequence from iterable containing {i!r}.")

        yield array[i]
