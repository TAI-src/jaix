from jaix.env.utils.switching_pattern.switching_pattern import (
    SeqRegSwitchingPatternConfig,
    SeqRegSwitchingPattern,
    SeqForcedSwitchingPatternConfig,
    SeqForcedSwitchingPattern,
)
from ttex.config import ConfigurableObjectFactory as COF


def test_init_reg():
    config = SeqRegSwitchingPatternConfig(wait_period=0.5)
    sp = COF.create(SeqRegSwitchingPattern, config, num_choices=3)

    assert sp.wait_period == 0.5
    assert sp.num_choices == 3
    assert hasattr(sp, "switch")


def test_switch_reg():
    config = SeqRegSwitchingPatternConfig(wait_period=0.5)
    sp = COF.create(SeqRegSwitchingPattern, config, num_choices=3)

    for _ in range(2):
        sp.reset()
        assert sp.switch(0.4) == 0
        assert sp.switch(0.5) == 1
        assert sp.switch(0.6) == 1
        assert sp.switch(1.3) == 2
        assert sp.switch(1.6) == -1


def test_switch_reg_force():
    config = SeqRegSwitchingPatternConfig(wait_period=0.5)
    sp = COF.create(SeqRegSwitchingPattern, config, num_choices=3)

    for _ in range(2):
        sp.reset()
        assert sp.switch(0.4) == 0
        assert sp.switch(0.4, [False, True, True]) == 1
        assert sp.switch(0.5) == 1
        assert sp.switch(0.6) == 1
        assert sp.switch(0.6, [False, False, True]) == 2
        assert sp.switch(0.99) == 2
        assert sp.switch(1.09) == 2
        assert sp.switch(1.1) == -1


def test_switch_reg_force_carry():
    config = SeqRegSwitchingPatternConfig(wait_period=0.5, carry_over=True)
    sp = COF.create(SeqRegSwitchingPattern, config, num_choices=3)

    for _ in range(2):
        sp.reset()
        assert sp.switch(0.4) == 0
        assert sp.switch(0.4, [False, True, True]) == 1
        assert sp.switch(0.5) == 1
        assert sp.switch(0.6) == 1
        assert sp.switch(0.6, [False, False, True]) == 2
        assert sp.switch(0.99, [False, False, True]) == 2
        assert sp.switch(0.99) == 1
        assert sp.switch(1.49, [False, False, True]) == 2
        assert sp.switch(1.5) == -1


def test_init_force():
    config = SeqForcedSwitchingPatternConfig()
    sp = COF.create(SeqForcedSwitchingPattern, config, num_choices=3)

    assert sp.num_choices == 3
    assert hasattr(sp, "switch")


def test_switch_force():
    config = SeqForcedSwitchingPatternConfig()
    sp = COF.create(SeqForcedSwitchingPattern, config, num_choices=3)

    for _ in range(2):
        sp.reset()
        assert sp.switch(0.4, [True, True, True]) == 0
        assert sp.switch(0.5, [True, True, True]) == 0
        assert sp.switch(1.6, [True, True, True]) == 0
        assert sp.switch(2, [False, True, True]) == 1
        assert sp.switch(2, [False, False, True]) == 2
        assert sp.switch(2, [False, False, True]) == 2
        assert sp.switch(2, [False, False, False]) == -1
