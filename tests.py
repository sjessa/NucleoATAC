#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import pytest
import sys

sys.exit(pytest.main(['-v', 'tests/']))
