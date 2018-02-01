#!/usr/bin/env python
#!/usr/bin/env pytest


from mean_shift import Mean_Shift


def test_Mean_Shift():
    a = Mean_Shift(3)
    x = [1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24]
    assert a.peak(x) == {0: 3, 1: 11, 2: 22}
