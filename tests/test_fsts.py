from dataclasses import dataclass

import fstlib
import numpy as np
import pytest

import medicc


@dataclass
class PairTest:
    seq_in: str
    seq_out: str
    expected_score_asym_wgd: float
    expected_score_asym_no_wgd: float
    expected_score_asym_total_cn: float
    expected_score_sym_wgd: float
    expected_score_sym_no_wgd: float
    expected_score_sym_total_cn: float


PAIR_TESTS = [
    PairTest('11111X1111', '22022X1111', 2, 2, 2, 2, 2, 2), # pair 0
    PairTest('11111X1111', '22022X2222', 2, 3, 3, 2, 3, 3), # pair 1
    PairTest('11101X1111', '10111X1111', np.inf, np.inf, np.inf, 2, 2, 2), # pair 2
    PairTest('33233X1111', '00000X1111', 3, 3, 3, 3, 3, 3), # pair 3
    PairTest('1111111111X1111111111', '2212222222X2222222222', 2, 3, 3, 2, 3, 3), # pair 4
    PairTest('2222222222X2222222222', '3323333333X3333323333', 3, 4, 4, 3, 4, 4), # pair 5
    PairTest('1111111X11X11X1111X1111', '3322112X22X23X2222X2200', 5, 9, 8, 5, 9, 8), # pair 6
    PairTest('1111111111X111X111', '3332222221X333X333', 4, 6, 3, 4, 6, 3), # pair 7
    PairTest('22222X222X222', '44444X444X444', 2, 6, 1, 2, 6, 1), # pair 8
    PairTest('2222222222', '2211001122', 2, 2, 2, 2, 2, 2), # pair 9
    PairTest('11X11X11X11', '22X33X33X33', 3, 7, 2, 3, 7, 2), # pair 10
    # PairTest('22X1X1X1', '10X2X2X2', 4, 5, 5, 4, 5, 5) # pair 11: Does not pass currently with WGD
]

FST_ASYMM_WGD = medicc.io.read_fst(no_wgd=False, total_copy_numbers=False, wgd_x2=False)
FST_ASYMM_NOWGD = medicc.io.read_fst(no_wgd=True, total_copy_numbers=False, wgd_x2=False)
FST_ASYMM_TOTAL = medicc.io.read_fst(no_wgd=False, total_copy_numbers=True, wgd_x2=False)

def _run_pair_test(pair_test: PairTest, is_wgd: bool, is_sym: bool, is_total_cn: bool) -> bool:
    """Runs individual pair test"""
    # maxcn = 8  # alphabet = 012345678; maxcn losses, maxcn-1 gains
    # sep = "X"
    # symbol_table = medicc.create_symbol_table(maxcn, sep)
    # fst = medicc.create_copynumber_fst(symbol_table, sep, enable_wgd=is_wgd)
    #fst = medicc.io.read_fst(no_wgd=(not is_wgd), total_copy_numbers=is_total_cn, wgd_x2=False)
    if is_wgd:
        if is_total_cn:
            fst = FST_ASYMM_TOTAL
        else:
            fst = FST_ASYMM_WGD
    else:
        fst = FST_ASYMM_NOWGD

    td = fstlib.factory.from_string(pair_test.seq_in, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    tg = fstlib.factory.from_string(pair_test.seq_out, isymbols=fst.input_symbols(), osymbols=fst.output_symbols(), arc_type=fst.arc_type())
    if is_sym:
        test_score = float(fstlib.kernel_score(fst, td, tg))
        if is_total_cn:
            expected_score = pair_test.expected_score_sym_total_cn
        elif is_wgd:
            expected_score = pair_test.expected_score_sym_wgd
        else:
            expected_score = pair_test.expected_score_sym_no_wgd
    else:
        test_score = float(fstlib.score(fst, td, tg))
        if is_total_cn:
            expected_score = pair_test.expected_score_asym_total_cn
        elif is_wgd:
            expected_score = pair_test.expected_score_asym_wgd
        else:
            expected_score = pair_test.expected_score_asym_no_wgd
    return test_score, expected_score


@pytest.mark.parametrize("pair", PAIR_TESTS)
def test_fstlib_sym_with_wgd(pair: PairTest):

    test_score, expected_score = _run_pair_test(pair, is_wgd=True, is_sym=True, is_total_cn=False)
    assert expected_score == test_score, f"expected: {expected_score}, test: {test_score}"

@pytest.mark.parametrize("pair", PAIR_TESTS)
def test_fstlib_sym_without_wgd(pair: PairTest):

    test_score, expected_score = _run_pair_test(pair, is_wgd=False, is_sym=True, is_total_cn=False)
    assert expected_score == test_score, f"expected: {expected_score}, test: {test_score}"

@pytest.mark.parametrize("pair", PAIR_TESTS)
def test_fstlib_asym_with_wgd(pair: PairTest):

    test_score, expected_score = _run_pair_test(pair, is_wgd=True, is_sym=False, is_total_cn=False)
    assert expected_score == test_score, f"expected: {expected_score}, test: {test_score}"

@pytest.mark.parametrize("pair", PAIR_TESTS)
def test_fstlib_asym_without_wgd(pair: PairTest):

    test_score, expected_score = _run_pair_test(pair, is_wgd=False, is_sym=False, is_total_cn=False)
    assert expected_score == test_score, f"expected: {expected_score}, test: {test_score}"

@pytest.mark.parametrize("pair", PAIR_TESTS)
def test_fstlib_asym_total_cn(pair: PairTest):

    test_score, expected_score = _run_pair_test(pair, is_wgd=True, is_sym=False, is_total_cn=True)
    assert expected_score == test_score, f"expected: {expected_score}, test: {test_score}"

@pytest.mark.parametrize("pair", PAIR_TESTS)
def test_fstlib_sym_total_cn(pair: PairTest):

    test_score, expected_score = _run_pair_test(pair, is_wgd=True, is_sym=True, is_total_cn=True)
    assert expected_score == test_score, f"expected: {expected_score}, test: {test_score}"
