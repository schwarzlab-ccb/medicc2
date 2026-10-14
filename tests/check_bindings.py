"""Exercise the refactored ABI boundary and old import compatibility."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from medicc.full_fst import ancestors, integration, core, pipeline
from medicc.full_fst.model import Model
from medicc.full_fst.pair_scorer import Model as PairScorer


def must_reject(action):
    try:
        action()
    except ValueError:
        return
    raise AssertionError('Expected ValueError')


def main():
    assert integration.FULL_TWO is ancestors.FULL_TWO
    assert integration.model is ancestors.model
    assert integration.reconstruct is ancestors.reconstruct
    assert pipeline.reconstruct is core.reconstruct
    must_reject(lambda: PairScorer(backend='missing'))
    for backend in ['direct_pair', 'reference_pair']:
        with PairScorer(backend=backend) as scorer:
            assert scorer.score('00014', '13207') == 5
            must_reject(lambda: scorer.score('1\0', '11'))
            must_reject(lambda: scorer.score('9', '1'))
            scorer.clear()
            assert scorer.stats() == [0] * 6
        scorer.close()
        must_reject(lambda: scorer.score('1', '2'))
    model = Model()
    left = model.string('1X1')
    right = model.string('8X8')
    try:
        retained = left.choose()
        right.choose()
        # Strings copied from native thread-local buffers survive subsequent calls.
        assert retained == ('1X1', 0)
        message = model.join(left, right)
        try:
            word, cost = model.resolve('1X1', message)
            assert isinstance(word, str) and isinstance(cost, int)
        finally:
            message.close()
        message.close()
        must_reject(message.choose)
    finally:
        left.close()
        right.close()
        model.close()
    model.close()
    must_reject(lambda: model.string('1'))
    print('PASS: shared compatibility imports, both pair ABIs, ownership, copied results, closed/error guards')


if __name__ == '__main__':
    main()
