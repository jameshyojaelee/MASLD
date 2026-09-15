#!/usr/bin/env python3
"""Matched from-scratch CNN and transformer regulatory controls."""

from __future__ import annotations

import random
from typing import Any, Mapping

import numpy as np


INPUT_LENGTH = 2114
OUTPUT_LENGTH = 1000
MODEL_IDS = ("sequence_cnn_control", "sequence_transformer_control")


class SequenceControlArchitectureError(ValueError):
    """Raised when a sequence-control architecture does not meet its frozen requirements."""


def _integer(
    parameters: Mapping[str, Any], key: str, *, minimum: int = 1
) -> int:
    try:
        value = int(parameters[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise SequenceControlArchitectureError(f"invalid integer parameter: {key}") from exc
    if value < minimum:
        raise SequenceControlArchitectureError(f"invalid integer parameter: {key}")
    return value


def _floating(
    parameters: Mapping[str, Any], key: str, *, minimum: float, maximum: float
) -> float:
    try:
        value = float(parameters[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise SequenceControlArchitectureError(f"invalid float parameter: {key}") from exc
    if not minimum <= value <= maximum:
        raise SequenceControlArchitectureError(f"invalid float parameter: {key}")
    return value


def _shared_parameters(parameters: Mapping[str, Any]) -> tuple[str, int, int, float]:
    model_id = str(parameters.get("control_model_id", ""))
    if model_id not in MODEL_IDS:
        raise SequenceControlArchitectureError("control_model_id is not admitted")
    input_length = _integer(parameters, "inputlen")
    output_length = _integer(parameters, "outputlen")
    if input_length != INPUT_LENGTH or output_length != OUTPUT_LENGTH:
        raise SequenceControlArchitectureError("sequence-control geometry differs")
    if (input_length - output_length) % 2:
        raise SequenceControlArchitectureError("profile crop is asymmetric")
    counts_loss_weight = _floating(
        parameters, "counts_loss_weight", minimum=0.0, maximum=1.0e9
    )
    return model_id, input_length, output_length, counts_loss_weight


def _seed_all(seed: int) -> None:
    import tensorflow as tf

    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)


def _build_cnn(parameters: Mapping[str, Any], input_length: int, output_length: int):
    import tensorflow as tf

    channels = _integer(parameters, "control_channels")
    blocks = _integer(parameters, "control_blocks")
    kernel_size = _integer(parameters, "control_kernel_size")
    stem_kernel_size = _integer(parameters, "control_stem_kernel_size")
    dropout = _floating(parameters, "control_dropout", minimum=0.0, maximum=0.9)
    if kernel_size % 2 != 1 or stem_kernel_size % 2 != 1:
        raise SequenceControlArchitectureError("CNN kernels must be odd")

    inputs = tf.keras.layers.Input(shape=(input_length, 4), name="sequence")
    values = tf.keras.layers.Conv1D(
        channels,
        stem_kernel_size,
        padding="same",
        activation="relu",
        name="cnn_stem",
    )(inputs)
    for index in range(blocks):
        residual = values
        values = tf.keras.layers.LayerNormalization(
            epsilon=1.0e-5, name=f"cnn_block_{index}_norm"
        )(values)
        values = tf.keras.layers.Conv1D(
            channels,
            kernel_size,
            padding="same",
            activation="relu",
            name=f"cnn_block_{index}_conv_1",
        )(values)
        values = tf.keras.layers.Dropout(
            dropout, name=f"cnn_block_{index}_dropout"
        )(values)
        values = tf.keras.layers.Conv1D(
            channels,
            kernel_size,
            padding="same",
            name=f"cnn_block_{index}_conv_2",
        )(values)
        values = tf.keras.layers.Add(name=f"cnn_block_{index}_residual")(
            [residual, values]
        )
        values = tf.keras.layers.Activation(
            "relu", name=f"cnn_block_{index}_activation"
        )(values)

    crop = (input_length - output_length) // 2
    profile_features = tf.keras.layers.Cropping1D(
        (crop, crop), name="cnn_profile_center_crop"
    )(values)
    profile = tf.keras.layers.Conv1D(
        1, 1, padding="same", name="cnn_profile_projection"
    )(profile_features)
    profile = tf.keras.layers.Flatten(name="logits_profile_predictions")(profile)
    pooled = tf.keras.layers.GlobalAveragePooling1D(name="cnn_global_average")(values)
    count = tf.keras.layers.Dense(1, name="logcount_predictions")(pooled)
    return tf.keras.Model(inputs=inputs, outputs=[profile, count], name="sequence_cnn_control")


def _build_transformer(
    parameters: Mapping[str, Any], input_length: int, output_length: int
):
    import tensorflow as tf

    width = _integer(parameters, "control_width")
    blocks = _integer(parameters, "control_blocks")
    heads = _integer(parameters, "control_heads")
    ffn_width = _integer(parameters, "control_ffn_width")
    token_stride = _integer(parameters, "control_token_stride")
    stem_kernel_size = _integer(parameters, "control_stem_kernel_size")
    profile_kernel_size = _integer(parameters, "control_profile_kernel_size")
    dropout = _floating(parameters, "control_dropout", minimum=0.0, maximum=0.9)
    attention_dropout = _floating(
        parameters, "control_attention_dropout", minimum=0.0, maximum=0.9
    )
    if width % heads or stem_kernel_size % 2 != 1 or profile_kernel_size % 2 != 1:
        raise SequenceControlArchitectureError("transformer width or kernels differ")
    if input_length % token_stride not in {0, 1, 2, 3}:
        raise SequenceControlArchitectureError("token stride is unsupported")

    inputs = tf.keras.layers.Input(shape=(input_length, 4), name="sequence")
    values = tf.keras.layers.Conv1D(
        width,
        stem_kernel_size,
        strides=token_stride,
        padding="same",
        activation="relu",
        name="transformer_tokenizer",
    )(inputs)
    values = tf.keras.layers.Conv1D(
        width,
        3,
        padding="same",
        name="transformer_relative_position_conv",
    )(values)
    for index in range(blocks):
        normalized = tf.keras.layers.LayerNormalization(
            epsilon=1.0e-5, name=f"transformer_block_{index}_attention_norm"
        )(values)
        attended = tf.keras.layers.MultiHeadAttention(
            num_heads=heads,
            key_dim=width // heads,
            dropout=attention_dropout,
            name=f"transformer_block_{index}_attention",
        )(normalized, normalized)
        attended = tf.keras.layers.Dropout(
            dropout, name=f"transformer_block_{index}_attention_dropout"
        )(attended)
        values = tf.keras.layers.Add(
            name=f"transformer_block_{index}_attention_residual"
        )([values, attended])
        normalized = tf.keras.layers.LayerNormalization(
            epsilon=1.0e-5, name=f"transformer_block_{index}_ffn_norm"
        )(values)
        feed_forward = tf.keras.layers.Dense(
            ffn_width, activation="relu", name=f"transformer_block_{index}_ffn_1"
        )(normalized)
        feed_forward = tf.keras.layers.Dropout(
            dropout, name=f"transformer_block_{index}_ffn_dropout"
        )(feed_forward)
        feed_forward = tf.keras.layers.Dense(
            width, name=f"transformer_block_{index}_ffn_2"
        )(feed_forward)
        values = tf.keras.layers.Add(name=f"transformer_block_{index}_ffn_residual")(
            [values, feed_forward]
        )

    count_features = tf.keras.layers.LayerNormalization(
        epsilon=1.0e-5, name="transformer_count_norm"
    )(values)
    pooled = tf.keras.layers.GlobalAveragePooling1D(
        name="transformer_global_average"
    )(count_features)
    count = tf.keras.layers.Dense(1, name="logcount_predictions")(pooled)

    profile_features = tf.keras.layers.UpSampling1D(
        size=token_stride, name="transformer_profile_upsample"
    )(values)
    upsampled_length = int(profile_features.shape[1])
    excess = upsampled_length - input_length
    if excess < 0 or excess % 2:
        raise SequenceControlArchitectureError("transformer upsampling geometry differs")
    if excess:
        profile_features = tf.keras.layers.Cropping1D(
            (excess // 2, excess // 2), name="transformer_upsample_crop"
        )(profile_features)
    crop = (input_length - output_length) // 2
    profile_features = tf.keras.layers.Cropping1D(
        (crop, crop), name="transformer_profile_center_crop"
    )(profile_features)
    profile = tf.keras.layers.Conv1D(
        1,
        profile_kernel_size,
        padding="same",
        name="transformer_profile_projection",
    )(profile_features)
    profile = tf.keras.layers.Flatten(name="logits_profile_predictions")(profile)
    return tf.keras.Model(
        inputs=inputs, outputs=[profile, count], name="sequence_transformer_control"
    )


def getModelGivenModelOptionsAndWeightInits(args, model_params):
    """Build the requested control for the native ChromBPNet trainer."""
    import tensorflow as tf

    from chrombpnet.training.utils.losses import multinomial_nll

    model_id, input_length, output_length, counts_loss_weight = _shared_parameters(
        model_params
    )
    seed = int(args.seed)
    if seed < 0:
        raise SequenceControlArchitectureError("seed must be nonnegative")
    _seed_all(seed)
    if model_id == "sequence_cnn_control":
        model = _build_cnn(model_params, input_length, output_length)
    else:
        model = _build_transformer(model_params, input_length, output_length)
    if model.input_shape != (None, INPUT_LENGTH, 4) or [
        tuple(shape) for shape in model.output_shape
    ] != [(None, OUTPUT_LENGTH), (None, 1)]:
        raise SequenceControlArchitectureError("built model geometry differs")
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=float(args.learning_rate)),
        loss=[multinomial_nll, "mse"],
        loss_weights=[1.0, counts_loss_weight],
    )
    return model


def save_model_without_bias(model, output_prefix) -> None:
    """The controls have no bias branch; the primary checkpoint is already bias-free."""
    del model, output_prefix

