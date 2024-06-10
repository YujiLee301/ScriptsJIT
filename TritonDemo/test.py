"""
Testing script is modified from 
https://github.com/triton-inference-server/server/blob/r21.02/qa/L0_backend_identity/identity_test.py
"""
import argparse
import numpy as np
import os
import re
import sys
import tritonclient.grpc as grpcclient
import tritonclient.http as httpclient
from tritonclient.utils import np_to_triton_dtype

FLAGS = None

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',
                        '--verbose',
                        action="store_true",
                        required=False,
                        default=False,
                        help='Enable verbose output')
    parser.add_argument('-u',
                        '--url',
                        type=str,
                        required=False,
                        help='Inference server URL.')
    parser.add_argument('-m',
                        '--model',
                        type=str,
                        required=False,
                        default='add_sub',
                        help='Model for inference')
    parser.add_argument(
        '-p',
        '--protocol',
        type=str,
        required=False,
        default='http',
        help='Protocol ("http"/"grpc") used to ' +
        'communicate with inference service. Default is "http".')

    FLAGS = parser.parse_args()
    if (FLAGS.protocol != "http") and (FLAGS.protocol != "grpc"):
        print("unexpected protocol \"{}\", expects \"http\" or \"grpc\"".format(
            FLAGS.protocol))
        exit(1)

    client_util = httpclient if FLAGS.protocol == "http" else grpcclient

    if FLAGS.url is None:
        FLAGS.url = "localhost:8000" if FLAGS.protocol == "http" else "localhost:8001"
    print("Flags url ", FLAGS.url)

    model_name = FLAGS.model
    assert model_name in ['particlenet_PT', 'add_sub'], "Error: model not in the repository."

    requests = 5
    with client_util.InferenceServerClient(FLAGS.url,
                                           verbose=FLAGS.verbose) as client:
        inputs = []
        results = []

        if model_name == 'add_sub':
            for i in range(requests):
                temp_0 = np.random.randn(16).astype(np.float32)
                temp_1 = np.random.randn(16).astype(np.float32)

                input_0 = client_util.InferInput('INPUT0', temp_0.shape, np_to_triton_dtype(np.float32))
                input_0.set_data_from_numpy(temp_0)

                input_1 = client_util.InferInput('INPUT1', temp_1.shape, np_to_triton_dtype(np.float32))
                input_1.set_data_from_numpy(temp_1)

                results.append(client.infer(model_name, [input_0, input_1]))

                inputs.append([temp_0, temp_1])

            for i in range(requests):
                print("\n*******\n Result: ")
                result = results[i]

                output_data0 = result.as_numpy("OUTPUT0")
                output_data1 = result.as_numpy("OUTPUT1")
                if output_data0 is None or output_data1 is None:
                    print("error: expected output missing")
                    sys.exit(1)

                print("input0: ", inputs[i][0])
                print("input1: ", inputs[i][1])
                print("sum: ", output_data0)
                print("sub: ", output_data1)

        elif model_name == 'particlenet_PT':
            for i in range(requests):
                temp_pf_points = np.random.randn(1,2,50).astype(np.float32)
                temp_pf_features = np.random.randn(1,25, 50).astype(np.float32)
                temp_pf_mask = np.random.randn(1,1,50).astype(np.float32)
                temp_sv_points = np.random.randn(1,2,10).astype(np.float32)
                temp_sv_features = np.random.randn(1,11,10).astype(np.float32)
                temp_sv_mask = np.random.randn(1,1,10).astype(np.float32)

                input_pf_points = client_util.InferInput("pf_points__0", temp_pf_points.shape, np_to_triton_dtype(np.float32))
                input_pf_points.set_data_from_numpy(temp_pf_points)
                input_pf_features = client_util.InferInput('pf_features__1', temp_pf_features.shape, np_to_triton_dtype(np.float32))
                input_pf_features.set_data_from_numpy(temp_pf_features)
                input_pf_mask = client_util.InferInput('pf_mask__2', temp_pf_mask.shape, np_to_triton_dtype(np.float32))
                input_pf_mask.set_data_from_numpy(temp_pf_mask)
                input_sv_points = client_util.InferInput('sv_points__3', temp_sv_points.shape, np_to_triton_dtype(np.float32))
                input_sv_points.set_data_from_numpy(temp_sv_points)
                input_sv_features = client_util.InferInput('sv_features__4', temp_sv_features.shape, np_to_triton_dtype(np.float32))
                input_sv_features.set_data_from_numpy(temp_sv_features)
                input_sv_mask = client_util.InferInput('sv_mask__5', temp_sv_mask.shape, np_to_triton_dtype(np.float32))
                input_sv_mask.set_data_from_numpy(temp_sv_mask)

                results.append(client.infer(model_name, [input_pf_points, input_pf_features, input_pf_mask, input_sv_points, input_sv_features, input_sv_mask]))
                inputs.append([temp_pf_points, temp_pf_features, temp_pf_mask, temp_sv_points, temp_sv_features, temp_sv_mask])

            for i in range(requests):
                print("\n*******\n Result: ")
                result = results[i]
                output_data0 = result.as_numpy("softmax__0")

                if output_data0 is None:
                    print("error: expected output missing")
                    sys.exit(1)

                print("output: ", output_data0)


    print("Passed all tests!")

