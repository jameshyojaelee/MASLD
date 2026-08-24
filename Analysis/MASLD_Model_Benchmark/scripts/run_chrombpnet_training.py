#!/usr/bin/env python3
"""Run native ChromBPNet training with all stochastic runtimes seeded first."""

from __future__ import annotations

import os
import random


def main() -> None:
    from chrombpnet.training.utils import argmanager

    arguments = argmanager.fetch_train_args()
    os.environ["PYTHONHASHSEED"] = str(arguments.seed)

    import numpy as np
    import tensorflow as tf

    random.seed(arguments.seed)
    np.random.seed(arguments.seed)
    tf.random.set_seed(arguments.seed)

    from chrombpnet.training import train

    train.main(arguments)


if __name__ == "__main__":
    main()
