import os
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")

import sys
import numpy as np
import pandas as pd


def _legacy_load_keras_h5_with_patch(h5_path):
    import h5py
    import tensorflow as tf
    from tensorflow.python.keras.models import model_from_json
    from tensorflow.python.keras import regularizers, initializers, constraints

    with h5py.File(h5_path, "r") as f:
        mc = f.attrs.get("model_config", None)
        if mc is None:
            if "model_config" in f:
                mc = f["model_config"][()]
            else:
                raise ValueError("模型文件中未找到 'model_config'，无法兼容性加载。")

        if isinstance(mc, bytes):
            mc = mc.decode("utf-8")
        elif not isinstance(mc, str):
            mc = str(mc)

        patched = mc.replace('"batch_shape"', '"batch_input_shape"')

    def _L2_wrapper(**kwargs):
        val = kwargs.get("l", kwargs.get("l2", 0.0))
        return regularizers.l2(val)

    def _L1_wrapper(**kwargs):
        val = kwargs.get("l", kwargs.get("l1", 0.0))
        return regularizers.l1(val)

    def _L1L2_wrapper(**kwargs):
        l1 = kwargs.get("l1", kwargs.get("l", 0.0))
        l2 = kwargs.get("l2", 0.0)
        return regularizers.l1_l2(l1=l1, l2=l2)

    custom_objects = {
        "L2": _L2_wrapper,
        "L1": _L1_wrapper,
        "L1L2": _L1L2_wrapper,
        "l2": _L2_wrapper,
        "l1": _L1_wrapper,
        "l1_l2": _L1L2_wrapper,

        "GlorotUniform": initializers.glorot_uniform(),
        "HeNormal": initializers.he_normal(),
        "MaxNorm": constraints.max_norm,
    }

    try:
        custom_objects.setdefault("gelu", getattr(tf.nn, "gelu", None))
    except Exception:
        pass

    def _swish(x):
        import tensorflow as tf
        return x * tf.nn.sigmoid(x)
    custom_objects.setdefault("swish", _swish)

    model = model_from_json(patched, custom_objects=custom_objects)

    try:
        model.load_weights(h5_path)
    except Exception:
        model.load_weights(h5_path, by_name=True, skip_mismatch=True)

    return model


def load_and_predict(model_file, X):
    ext = os.path.splitext(str(model_file))[1].lower()

    if ext in {".joblib", ".pkl"}:
        import joblib
        model = joblib.load(model_file)
        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(X)
            return proba[:, 1]
        elif hasattr(model, "decision_function"):
            z = model.decision_function(X)
            return 1.0 / (1.0 + np.exp(-z))
        else:
            y = model.predict(X)
            return np.asarray(y, dtype=float).ravel()

    elif ext == ".h5":
        model = _legacy_load_keras_h5_with_patch(model_file)
        y_pred = model.predict(X, verbose=0)
        y_pred = np.asarray(y_pred)
        if y_pred.ndim == 2 and y_pred.shape[1] == 1:
            y_proba = y_pred[:, 0]
        elif y_pred.ndim == 2 and y_pred.shape[1] >= 2:
            y_proba = y_pred[:, -1]
        else:
            y_proba = y_pred.squeeze()
        return y_proba.ravel()

    elif ext == ".keras":
        raise RuntimeError(
            ".keras 单文件不受此旧版 TensorFlow/Keras 支持；"
            "请改用已转换的 .h5，或在较新环境加载。"
        )
    else:
        raise ValueError("不支持的模型文件类型：{}（仅支持 .joblib/.pkl/.h5/.keras）".format(ext))


def predict_with_model(input_file, model_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    X = df.iloc[:, 1:].values  

    y_proba = load_and_predict(model_file, X)
    labels = np.where(y_proba >= 0.5, "Virus", "Not Virus")

    out = pd.DataFrame({
        "ID": df.iloc[:, 0].values,
        "Score": y_proba,
        "Prediction": labels
    })
    out.to_csv(output_file, sep='\t', index=False)
    print("✅ Prediction completed. Results saved to:", output_file)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python predict.py <input_file> <model_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    model_file = sys.argv[2]
    output_file = sys.argv[3]
    predict_with_model(input_file, model_file, output_file)
