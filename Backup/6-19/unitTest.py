#!/usr/bin/python -tt
# Copyright Nima Emami, 2012

import unittest

class TestSequenceFunctions(unittest.TestCase):

    def test_choice(self):
        element = random.choice(self.seq)
        self.assertTrue(element in self.seq)
        # should raise an exception for an immutable sequence
        self.assertRaises(TypeError, random.shuffle, (1,2,3))

    def test_sample(self):
        with self.assertRaises(ValueError):
            random.sample(self.seq, 20)
        for element in random.sample(self.seq, 5):
            self.assertTrue(element in self.seq)

if __name__ == '__main__':
    unittest.main()