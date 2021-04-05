using System;
using System.Collections.Generic;
using Numpy;
using System.Linq;

namespace ConsoleGpsConv {

    public static class Simil {

        public static NDarray _get_scalar(NDarray alpha_0, NDarray q_coords = null) {
            NDarray scalar;
            //var temp = new NDarray[] { alpha_0 };
            if (q_coords == null) {
                //scalar = np.einsum("i->", temp);
                scalar = np.sum(alpha_0, axis: -1);
            } else {
                //scalar = np.einsum("i,i->", temp, (q_coords * q_coords).sum(0));
                scalar = np.dot(alpha_0, (q_coords * q_coords).sum(0));
            }
            return scalar;
        }

        public static NDarray _get_q_matrix(NDarray quaternions) {
            var list = new List<float[]>();
            if (quaternions.ndim == 1)
            {
                var floats = new List<float>();
                for (int i = 0; i < quaternions.len; i++)
                {
                    floats.Add(getFloatValue(quaternions[i]));
                }
                list.Add(floats.ToArray());
            }
            else
            {

                for (int i = 0; i < quaternions.len; i++)
                {
                    list.Add((quaternions[i]).GetData<float>());
                }

            }
            var q_matrix = (from q in list
                            select new float[][] {
                    new float[] {
                        q[3],
                        -q[2],
                        q[1],
                        q[0]
                    },
                    new float[] {
                        q[2],
                        q[3],
                        -q[0],
                        q[1]
                    },
                    new float[] {
                        -q[1],
                        q[0],
                        q[3],
                        q[2]
                    },
                    new float[] {
                        -q[0],
                        -q[1],
                        -q[2],
                        q[3]
                    }
                }).ToArray();

            var matrix = getMatrix(q_matrix);

            return np.array(matrix);
        }

        public static NDarray _get_w_matrix(NDarray quaternions) {
            var list = new List<float[]>();
            if(quaternions.ndim == 1)
            {
                var floats = new List<float>();
                for (int i = 0; i < quaternions.len; i++)
                {
                    floats.Add(getFloatValue(quaternions[i]));
                }
                list.Add(floats.ToArray());
            }
            else
            {

                for (int i = 0; i < quaternions.len; i++)
                {
                    list.Add((quaternions[i]).GetData<float>());
                }

            }
            var w_matrix = (from q in list
                            select new float[][] {
                    new float[] {
                        q[3],
                        q[2],
                        -q[1],
                        q[0]
                    },
                    new float[] {
                        -q[2],
                        q[3],
                        q[0],
                        q[1]
                    },
                    new float[] {
                        q[1],
                        -q[0],
                        q[3],
                        q[2]
                    },
                    new float[] {
                        -q[0],
                        -q[1],
                        -q[2],
                        q[3]
                    }
                }).ToArray();


            var matrix = getMatrix(w_matrix);

            return np.array(matrix);
        }

        public static float[,,] getMatrix(float[][][] matrix) {
            float[,,] test = new float[matrix.Count(), matrix.First().Length, 4];
            for (int i = 0; i < matrix.Count(); i++) {
                for (int j = 0; j < matrix[i].Length; j++) {
                    for (int k = 0; k < matrix[i][j].Length; k++) {
                        test[i, j, k] = matrix[i][j][k];
                    }
                }
            }
            return test;
        }

        public static NDarray _get_abc_matrices(NDarray alpha_0, NDarray m1, NDarray m2 = null) {
            NDarray matrix;
            if (m2 == null) {
                var temp = new NDarray[] { alpha_0, m1 };
                //matrix = np.einsum("i,ijk->jk", temp);
                matrix = np.tensordot(alpha_0, m1, new[] { 0, 0 });
            }
            else {
                var temp = new NDarray[] { alpha_0, np.matmul(np.transpose(m1, new[] { 0, 2, 1 }), m2) };
                //matrix = np.einsum("i,ijk->jk", temp);
                matrix = np.tensordot(alpha_0, np.matmul(np.transpose(m1, new[] { 0, 2, 1 }), m2), new[] { 0, 0 });
                    //np.tensordot(t,I, axes=([0],[0]))

            }
            return matrix;
        }

        public static NDarray _get_blc_matrix(NDarray b_matrix, double lambda_i, NDarray c_matrix) {
            var blc_matrix = b_matrix - lambda_i * c_matrix;
            return blc_matrix;
        }

        public static NDarray _get_d_matrix(double li, NDarray cs, NDarray am, NDarray blcm) {
            var d_matrix = 2 * li * am + (1 / cs) * np.matmul(blcm.T, blcm);
            return d_matrix;
        }


        public static (NDarray, NDarray) _get_r_quat(NDarray d_matrix) {
            NDarray eigvals, eigvects;
            (eigvals, eigvects) = np.linalg.eig(d_matrix);
            NDarray beta_1 = np.argmax(eigvals);
            //NDarray r_quat = eigvects[":, beta_1"];
            var r_quat = eigvects[":", beta_1];
            return (beta_1, r_quat);
        }

        public static float _get_lambda_next(NDarray am, NDarray bs, NDarray bm, NDarray cs, NDarray cm, NDarray rq)
        {
            var expr_1 = np.matmul(np.matmul(rq.T, am), rq);
            var expr_2 = (1 / cs) * np.matmul(np.matmul(np.matmul(rq.T, bm.T), cm), rq);
            var expr_3 = (1 / cs) * np.matmul(np.matmul(np.matmul(rq.T, cm.T), cm), rq);
            var lambda_next = (expr_1 - expr_2) / (bs - expr_3);
            return getFloatValue(lambda_next);
            //return lambda_next.astype(np.float64).GetData<float>()[0];
            //return float.Parse(lambda_next.astype(np.float64));

        }

        private static float getFloatValue(NDarray array)
        {
            float res = 0;
            if(!float.TryParse(array.real.repr, System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out res))
            {
                throw new Exception();
            }
            return res;
        }

        //def _get_lambda_next(am, bs, bm, cs, cm, rq):
        //expr_1 = rq.T @ am @ rq
        //expr_2 = (1 / cs) * (rq.T @ bm.T @ cm @ rq)
        //expr_3 = (1/cs) * (rq.T @ cm.T @ cm @ rq)
        //lambda_next = (expr_1-expr_2) / (bs-expr_3)
        //return lambda_next


        public static (NDarray, NDarray, NDarray, NDarray, float, float) _get_solution(NDarray am, NDarray bs, NDarray bm, NDarray cs, NDarray cm, bool scale, float li, float i) {
            var blc_matrix = _get_blc_matrix(bm, li, cm);
            var d_matrix = _get_d_matrix(li, cs, am, blc_matrix);
            NDarray beta_1, r_quat;
            (beta_1, r_quat) = _get_r_quat(d_matrix);
            if (scale == false) {
                return (blc_matrix, d_matrix, beta_1, r_quat, li, i);
            } else {
                var lambda_next = _get_lambda_next(am, bs, bm, cs, cm, r_quat);
                if (Math.Abs(li - lambda_next) < 0.000001)
                {
                    return (blc_matrix, d_matrix, beta_1, r_quat, li, i);
                }
                else
                {
                    li = lambda_next;
                    i = i + 1;
                    //li, i = lambda_next, i + 1;
                    return _get_solution(am, bs, bm, cs, cm, scale, li, i);
                }

                //lambda_next = _get_lambda_next(am, bs, bm, cs, cm, r_quat)
                //if abs(li - lambda_next) < 0.000001:
                //    return blc_matrix, d_matrix, beta_1, r_quat, li, i
                //else:
                //    li, i = lambda_next, i + 1
                //    return _get_solution(am, bs, bm, cs, cm, scale, li, i)
                //throw new Exception("not implemeted!");
            }
        }

        public static NDarray _get_r_matrix(NDarray r_quat) {
            //var temp = new NDarray[] { r_quat };
            var r_w_matrix = _get_w_matrix(r_quat)[0];
            var r_q_matrix = _get_q_matrix(r_quat)[0];
            var r_matrix = np.matmul(r_w_matrix.T, r_q_matrix)[":3,:3"];
            return r_matrix;
        }

        public static NDarray _get_s_quat(NDarray c_scalar, NDarray blcm, NDarray r_quat) {
            var s_quat = 1 / (2 * c_scalar) * np.matmul(blcm, r_quat);
            return s_quat;
        }

        public static NDarray _get_t_vector(NDarray r_quat, NDarray s_quat) {
            //var temp = new NDarray[] { r_quat };
            var r_w_matrix = _get_w_matrix(r_quat)[0];
            var t_vector = 2 * np.matmul(r_w_matrix.T, s_quat)[":3"];
            return t_vector;
        }



        //# ================
        //# Process function
        //# ================

        public static (float, NDarray, NDarray) process(float[,] source_points,
                    float[,] target_points,
                    NDarray alpha_0 = null,
                    bool scale = true,
                    float lambda_0 = 1.0f) {

            //"""
            //Find similarity transformation parameters given a set of control points


            //Parameters
            //----------
            //source_points : array_like
            //    The function will try to cast it to a numpy array with shape:
            //    ``(n, 3)``, where ``n`` is the number of points.
            //    Two points is the minimum requeriment (in that case, the solution
            //    will map well all points that belong in the rect that passes
            //    through both control points).
            //target_points : array_like
            //    The function will try to cast it to a numpy array with shape:
            //    ``(n, 3)``, where ``n`` is the number of points.
            //    The function will check that there are as many target points
            //    as source points.
            //alpha_0 : array_like, optional
            //    Per point weights.
            //    If provided, the function will try to cast to a numpy array with
            //    shape: ``(n,)``.
            //scale : boolean, optional
            //    Allow to find a multiplier factor different from lambda_0.
            //    Default is True.
            //lambda_0 : float, optional
            //    Multiplier factor to find the first solution. Default is 1.0.
            //    If `scale= True`, a recursion is implemented to find a better
            //    value. If it is negative, forces mirroring. Can't be zero.
            //Returns
            //-------
            //lambda_i : float
            //    Multiplier factor.
            //r_matrix : numpy.ndarray
            //    Rotation matrix.
            //t_vector : numpy.ndarray
            //    Translation (column) vector.
            //"""


            //declarations and checkups

            var source_coords = np.array(source_points, dtype: np.float64).T;

            if (source_coords.ndim != 2)
                throw new Exception("source_points array must have dimension = 2.");


            if (source_coords.shape[0] != 3)
                throw new Exception("There are not three coordinates in source points.");

            var n = source_coords.shape[1];

            //if (n == 1 || (source_coords[null, 0] == source_coords).all())
            //    throw new Exception("There are not two distinct source points.");


            var target_coords = np.array(target_points, dtype: np.float64).T;


            if (target_coords.ndim != 2)
                throw new Exception("target_points array must have dimension = 2.");

            if (target_coords.shape[0] != 3)
                throw new Exception("There are not three coordinates in target points.");

            if (target_coords.shape[1] != n)
                throw new Exception("There are not as many target points as source points.");

            if (alpha_0 == null)
                alpha_0 = np.ones(n);
            else
                alpha_0 = np.array(alpha_0, dtype: np.float64);

            if (alpha_0.ndim != 1)
                throw new Exception("alpha_0 array must have dimension = 1.");

            //if (alpha_0.shape != (n,))
            //    throw new Exception("There are not as many alpha_0 coefficients as control points.");

            //lambda_0 = float(lambda_0);


            if (lambda_0 == 0)
                throw new Exception("lambda_0 cannot be zero.");


            //processes

            var source_q_coords = np.concatenate((source_coords, np.zeros((1, n))));


            var target_q_coords = np.concatenate((target_coords, np.zeros((1, n))));

            var b_scalar = _get_scalar(alpha_0, source_q_coords);

            //var temp = new NDarray[] { alpha_0 };
            //var c_scalar = np.einsum("i->", temp);
            var c_scalar = np.sum(alpha_0, axis: -1);


            var q0_w_matrix = _get_w_matrix(source_q_coords.T);


            var qt_q_matrix = _get_q_matrix(target_q_coords.T);


            var a_matrix = _get_abc_matrices(alpha_0, q0_w_matrix, qt_q_matrix);


            var b_matrix = _get_abc_matrices(alpha_0, qt_q_matrix);


            var c_matrix = _get_abc_matrices(alpha_0, q0_w_matrix);

            float lambda_i = lambda_0;
            float i = 1;

            //var lambda_i, i = lambda_0, 1;

            NDarray blc_matrix, d_matrix, beta_1, r_quat;
            (blc_matrix, d_matrix, beta_1, r_quat, lambda_i, i) = _get_solution(a_matrix,
                                                                              b_scalar,
                                                                              b_matrix,
                                                                              c_scalar,
                                                                              c_matrix,
                                                                              scale,
                                                                              lambda_i,
                                                                              i);


            var r_matrix = _get_r_matrix(r_quat);


            var s_quat = _get_s_quat(c_scalar, blc_matrix, r_quat);


            var t_vector = np.array(_get_t_vector(r_quat, s_quat)).reshape(3, 1);


            return (lambda_i, r_matrix, t_vector);

        }

    }
}