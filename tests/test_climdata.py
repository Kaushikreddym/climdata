#!/usr/bin/env python

import unittest
import climdata

class TestClimdata(unittest.TestCase):
    def test_version(self):
        self.assertIsNotNone(climdata.__version__)

