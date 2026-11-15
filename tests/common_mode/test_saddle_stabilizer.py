"""Freezes tests/probes/saddle_stabilizer.py: on the empty field's saddle every row reads one cone."""


def _silent(*args, **kwargs):
    pass


def test_saddle_stabilizer_every_row_reads_one_cone():
    """The stabilizer of the (1,3) saddle has dimension 6 and leaves no
    vector, one rank-2 tensor (the saddle), no rank-3 tensor and four of
    rank 4.  A scalar row's symbol is a multiple of the inverse saddle.
    A vector row's determinant vanishes on its cone and not off it.  A
    Weyl row has one invariant first order symbol, whose determinant is
    a multiple of the same quadratic form (asserted inside the probe).
    Waves of the scalar, Weyl and vector rows weigh their energy and
    run at the saddle's speed."""
    from probes import saddle_stabilizer as ss

    r = ss.run(report=_silent)
    assert r["stabilizer_dim"] == 6
    assert r["invariant_counts"] == {"rank 1": 0, "rank 2": 1, "rank 3": 0, "rank 4": 4}
    assert r["vector_det"]["null"] < 1e-12 and r["vector_det"]["generic"] > 1e-3
    assert r["weyl_invariants"] == 1
    assert r["weight_over_energy"] < 1e-12 and r["weyl_vector_error"] < 1e-12


def test_saddle_stabilizer_a_frame_gives_each_row_its_own_speed():
    """The statement can fail.  With one vector in the vacuum the
    stabilizer falls to dimension 3, and a second symmetric rank-2
    tensor and a second Weyl symbol appear.  Scalar rows built from the
    two tensors run at speeds that differ from each other and from the
    saddle's.  Their weight falls below their energy by a different
    factor for each."""
    from probes import saddle_stabilizer as ss

    r = ss.run(report=_silent)
    assert r["frame_stabilizer_dim"] == 3 and r["frame_vectors"] == 1
    assert r["frame_symmetric_rank2"] == 2 and r["frame_weyl_invariants"] == 2
    speeds = r["frame_speeds"]
    assert max(speeds) - min(speeds) > 1e-2
    assert all(abs(s - r["shared_speed"]) > 1e-3 for s in speeds)
    assert max(r["frame_weight"]) < 1.0 and r["frame_eotvos"] > 1e-2
    assert max(r["frame_weight_weyl"] + r["frame_weight_vector"]) < 0.99
