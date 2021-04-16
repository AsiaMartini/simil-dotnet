using MathNet.Numerics.LinearAlgebra.Double;
using NumpyDotNet;
using System;
using System.Collections.Generic;
using System.Linq;

namespace ConsoleGpsConv {

    public static class Simil {

        public static ndarray _get_scalar(ndarray alpha_0, ndarray q_coords = null) {
            ndarray scalar;
            //var temp = new ndarray[] { alpha_0 };
            if (q_coords == null) {
                //scalar = np.einsum("i->", temp);
                scalar = np.sum(alpha_0, axis: -1);
            } else {
                //scalar = np.einsum("i,i->", temp, (q_coords * q_coords).sum(0));
                scalar = np.dot(alpha_0, (q_coords * q_coords).Sum(0));
            }
            return scalar;
        }

        public static ndarray _get_q_matrix(ndarray quaternions) {


            if (quaternions.ndim == 0) {

                var temp2 = new List<List<double>>();
                temp2.Add(new List<double> {
                        (double)quaternions[3],
                        -(double)quaternions[2],
                        (double)quaternions[1],
                        (double)quaternions[0]
                    });
                temp2.Add(new List<double> {
                        (double)quaternions[2],
                        (double)quaternions[3],
                        -(double)quaternions[0],
                        (double)quaternions[1]
                    });
                temp2.Add(new List<double> {
                        -(double)quaternions[1],
                        (double)quaternions[0],
                        (double)quaternions[3],
                        (double)quaternions[2]
                    });
                temp2.Add(new List<double> {
                         -(double)quaternions[0],
                        -(double)quaternions[1],
                        -(double)quaternions[2],
                        (double)quaternions[3]
                    });


                double[][] arrays = temp2.Select(a => a.ToArray()).ToArray();

                var matrix = getMatrix2D(arrays);

                return np.array(matrix);


            } else {

                List<List<List<double>>> result = new List<List<List<double>>>();
                for (int i = 0; i <= quaternions.ndim; i++) {
                    var q = (ndarray)quaternions[i];

                    var temp2 = new List<List<double>>();
                    temp2.Add(new List<double> {
                            (double)q[3],
                            -(double)q[2],
                            (double)q[1],
                            (double)q[0]
                        });
                    temp2.Add(new List<double> {
                            (double)q[2],
                            (double)q[3],
                            -(double)q[0],
                            (double)q[1]
                        });
                    temp2.Add(new List<double> {
                            -(double)q[1],
                            (double)q[0],
                            (double)q[3],
                            (double)q[2]
                        });
                    temp2.Add(new List<double> {
                            -(double)q[0],
                            -(double)q[1],
                            -(double)q[2],
                            (double)q[3]
                        });

                    result.Add(temp2);

                }

                double[][][] arrays = result.Select(a => a.Select(b => b.ToArray()).ToArray()).ToArray();

                var matrix = getMatrix(arrays);

                return np.array(matrix);
            }
        }

        public static ndarray _get_w_matrix(ndarray quaternions) {


            if (quaternions.ndim == 0) {

                var temp2 = new List<List<double>>();
                temp2.Add(new List<double> {
                        (double)quaternions[3],
                        (double)quaternions[2],
                        -(double)quaternions[1],
                        (double)quaternions[0]
                    });
                temp2.Add(new List<double> {
                        -(double)quaternions[2],
                        (double)quaternions[3],
                        (double)quaternions[0],
                        (double)quaternions[1]
                    });
                temp2.Add(new List<double> {
                        (double)quaternions[1],
                        -(double)quaternions[0],
                        (double)quaternions[3],
                        (double)quaternions[2]
                    });
                temp2.Add(new List<double> {
                         -(double)quaternions[0],
                        -(double)quaternions[1],
                        -(double)quaternions[2],
                        (double)quaternions[3]
                    });


                double[][] arrays = temp2.Select(a => a.ToArray()).ToArray();

                var matrix = getMatrix2D(arrays);

                return np.array(matrix);


            } else {


                List<List<List<double>>> result = new List<List<List<double>>>();
                for (int i = 0; i <= quaternions.ndim; i++) {
                    var q = (ndarray)quaternions[i];

                    var temp2 = new List<List<double>>();
                    temp2.Add(new List<double> {
                        (double)q[3],
                        (double)q[2],
                        -(double)q[1],
                        (double)q[0]
                    });
                    temp2.Add(new List<double> {
                        -(double)q[2],
                        (double)q[3],
                        (double)q[0],
                        (double)q[1]
                    });
                    temp2.Add(new List<double> {
                        (double)q[1],
                        -(double)q[0],
                        (double)q[3],
                        (double)q[2]
                    });
                    temp2.Add(new List<double> {
                         -(double)q[0],
                        -(double)q[1],
                        -(double)q[2],
                        (double)q[3]
                    });

                    result.Add(temp2);

                }

                double[][][] arrays = result.Select(a => a.Select(b => b.ToArray()).ToArray()).ToArray();

                var matrix = getMatrix(arrays);

                return np.array(matrix);
            }


        }

        public static double[,,] getMatrix(double[][][] matrix) {
            double[,,] test = new double[matrix.Count(), matrix.First().Length, 4];
            for (int i = 0; i < matrix.Count(); i++) {
                for (int j = 0; j < matrix[i].Length; j++) {
                    for (int k = 0; k < matrix[i][j].Length; k++) {
                        test[i, j, k] = matrix[i][j][k];
                    }
                }
            }
            return test;
        }

        public static double[,] getMatrix2D(double[][] matrix) {
            double[,] test = new double[matrix.Count(), matrix.First().Length];
            for (int i = 0; i < matrix.Count(); i++) {
                for (int j = 0; j < matrix[i].Length; j++) {
                    test[i, j] = matrix[i][j];
                }
            }
            return test;
        }


        public static ndarray _get_abc_matrices(ndarray alpha_0, ndarray m1, ndarray m2 = null) {
            ndarray matrix;
            if (m2 == null) {
                //matrix = np.einsum("i,ijk->jk", temp);
                matrix = np.tensordot(alpha_0, m1, axes: (new long[] { 0 }, new long[] { 0 }));
            }
            else {
                //matrix = np.einsum("i,ijk->jk", temp);
                var trans = np.transpose(m1, new long[] { 0, 2, 1 });
                var temp = np.matmul(trans, m2);
                //matrix = np.tensordot(alpha_0, temp);

                matrix = np.tensordot(alpha_0, temp, axes: (new long[] { 0 }, new long[] { 0 }));
            }
            return matrix;
        }

        public static ndarray _get_blc_matrix(ndarray b_matrix, double lambda_i, ndarray c_matrix) {
            var blc_matrix = b_matrix - lambda_i * c_matrix;
            return blc_matrix;
        }

        public static ndarray _get_d_matrix(double li, ndarray cs, ndarray am, ndarray blcm) {
            var d_matrix = 2 * li * am + (1 / cs) * np.matmul(blcm.T, blcm);
            return d_matrix;
        }


        public static (ndarray, ndarray) _get_r_quat2(ndarray d_matrix) {
            ndarray eigvals = null, eigvects = null;
            //(eigvals, eigvects) = np.linalg.eig(d_matrix);
            ndarray beta_1 = np.argmax(eigvals);
            var r_quat = eigvects[":", beta_1];
            return (beta_1, (ndarray)r_quat);
        }



        public static (ndarray, ndarray) _get_r_quat(ndarray d_matrix) {
            ndarray eigvals, eigvects;

            //original function
            //(eigvals, eigvects) = np.linalg.eig(d_matrix);

            List<List<double>> result = new List<List<double>>();
            for (int i = 0; i < d_matrix.shape.iDims[1]; i++) {
                var q = (ndarray)d_matrix[i];
                var temp2 = new List<double>() {
                        (double)q[0],
                        (double)q[1],
                        (double)q[2],
                        (double)q[3]
                    };
                result.Add(temp2);
            }


            double[][] arrays = result.Select(a => a.ToArray()).ToArray();

            var matrix = DenseMatrix.OfArray(getMatrix2D(arrays));
            var evd = matrix.Evd();
            Console.WriteLine("-----------------");
            Console.WriteLine($"\nevd.EigenValues :\n{evd.EigenValues}");
            Console.WriteLine($"\nevd.EigenVectors :\n{evd.EigenVectors}");

            var d = evd.EigenValues.Select(x => x.Real).ToArray<double>();
            eigvals = np.array(d);
            var d2 = evd.EigenVectors.ToArray();
            eigvects = np.array(d2);

            ndarray beta_1 = np.argmax(eigvals);
            Console.WriteLine($"\nbeta_1 :\n{beta_1}");
            var r_quat = eigvects[":", beta_1];
            return (beta_1, (ndarray)r_quat);

        }

        public static double _get_lambda_next(ndarray am, ndarray bs, ndarray bm, ndarray cs, ndarray cm, ndarray rq)
        {
            var expr_1 = np.matmul(np.matmul(rq.T, am), rq);
            var expr_2 = (1 / cs) * np.matmul(np.matmul(np.matmul(rq.T, bm.T), cm), rq);
            var expr_3 = (1 / cs) * np.matmul(np.matmul(np.matmul(rq.T, cm.T), cm), rq);
            var lambda_next = (expr_1 - expr_2) / (bs - expr_3);
            return getFloatValue((ndarray)lambda_next);
        }

        private static double getFloatValue(ndarray array)
        {
            double res = 0;
            if(!double.TryParse(array.Real.ToString(), System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture, out res))
            {
                throw new Exception();
            }
            return res;
        }


        public static (ndarray, ndarray, ndarray, ndarray, double, double) _get_solution(ndarray am, ndarray bs, ndarray bm, ndarray cs, ndarray cm, bool scale, double li, double i) {
            var blc_matrix = _get_blc_matrix(bm, li, cm);
            var d_matrix = _get_d_matrix(li, cs, am, blc_matrix);
            ndarray beta_1, r_quat;
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
                    return _get_solution(am, bs, bm, cs, cm, scale, li, i);
                }
            }
        }

        public static ndarray _get_r_matrix(ndarray r_quat) {
            var r_w_matrix = _get_w_matrix(r_quat);
            var r_q_matrix = _get_q_matrix(r_quat);
            var r_matrix = np.matmul(r_w_matrix.T, r_q_matrix)[":3",":3"];
            return (ndarray)r_matrix;
        }

        public static ndarray _get_s_quat(ndarray c_scalar, ndarray blcm, ndarray r_quat) {
            var s_quat = 1 / (2 * c_scalar) * np.matmul(blcm, r_quat);
            return s_quat;
        }

        public static ndarray _get_t_vector(ndarray r_quat, ndarray s_quat) {
            var r_w_matrix = _get_w_matrix(r_quat);
            var t_vector = 2 * (ndarray)np.matmul(r_w_matrix.T, s_quat)[":3"];
            return t_vector;
        }



        //# ================
        //# Process function
        //# ================

        public static (double, ndarray, ndarray) process(double[,] source_points,
                    double[,] target_points,
                    ndarray alpha_0 = null,
                    bool scale = true,
                    double lambda_0 = 1.0f) {

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
            //lambda_0 : double, optional
            //    Multiplier factor to find the first solution. Default is 1.0.
            //    If `scale= True`, a recursion is implemented to find a better
            //    value. If it is negative, forces mirroring. Can't be zero.
            //Returns
            //-------
            //lambda_i : double
            //    Multiplier factor.
            //r_matrix : numpy.ndarray
            //    Rotation matrix.
            //t_vector : numpy.ndarray
            //    Translation (column) vector.
            //"""


            //declarations and checkups

            var source_coords = np.array(source_points, dtype: np.Float64).T;

            if (source_coords.ndim != 2)
                throw new Exception("source_points array must have dimension = 2.");


            if (source_coords.shape.iDims[0] != 3)
                throw new Exception("There are not three coordinates in source points.");

            var n = source_coords.shape.iDims[1];

            //if (n == 1 || (source_coords[null, 0] == source_coords).all())
            //    throw new Exception("There are not two distinct source points.");


            var target_coords = np.array(target_points, dtype: np.Float64).T;


            if (target_coords.ndim != 2)
                throw new Exception("target_points array must have dimension = 2.");

            if (target_coords.shape.iDims[0] != 3)
                throw new Exception("There are not three coordinates in target points.");

            if (target_coords.shape.iDims[1] != n)
                throw new Exception("There are not as many target points as source points.");

            if (alpha_0 == null)
                alpha_0 = np.ones(n);
            else
                alpha_0 = np.array(alpha_0, dtype: np.Float64);

            if (alpha_0.ndim != 1)
                throw new Exception("alpha_0 array must have dimension = 1.");

            //if (alpha_0.shape != (n,))
            //    throw new Exception("There are not as many alpha_0 coefficients as control points.");

            //lambda_0 = double(lambda_0);


            if (lambda_0 == 0)
                throw new Exception("lambda_0 cannot be zero.");


            //processes

            var source_q_coords = np.concatenate((source_coords, np.zeros(new shape(1, n), np.Int32)));


            var target_q_coords = np.concatenate((target_coords, np.zeros(new shape(1, n), np.Int32)));

            var b_scalar = _get_scalar(alpha_0, source_q_coords);

            //var temp = new ndarray[] { alpha_0 };
            //var c_scalar = np.einsum("i->", temp);
            var c_scalar = np.sum(alpha_0, axis: -1);


            var q0_w_matrix = _get_w_matrix(source_q_coords.T);


            var qt_q_matrix = _get_q_matrix(target_q_coords.T);


            var a_matrix = _get_abc_matrices(alpha_0, q0_w_matrix, qt_q_matrix);


            var b_matrix = _get_abc_matrices(alpha_0, qt_q_matrix);


            var c_matrix = _get_abc_matrices(alpha_0, q0_w_matrix);

            double lambda_i = lambda_0;
            double i = 1;

            //var lambda_i, i = lambda_0, 1;

            ndarray blc_matrix, d_matrix, beta_1, r_quat;
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







        public static void test() {

            double[,] source_points = new double[,] {
                {0, 0, 0},
                {0, 2, 2},
                {2, 3, 1}
            };

            double[,] target_points = new double[,] {
                {3, 7, 5},
                {6, 7, 2},
                {4.5f, 4, 0.5f}
            };


            var source_coords = np.array(source_points, dtype: np.Float64).T;
            var target_coords = np.array(target_points, dtype: np.Float64).T;
            var n = source_coords.shape.iDims[1];
            var alpha_0 = np.ones(n);
            float lambda_0 = 1.0f;
            bool scale = true;


            //eq 29
            //var source_q_coords = np.concatenate((source_coords, np.zeros((1, n))));

            var source_q_coords = np.concatenate((source_coords, np.zeros(new shape(1, n), np.Float64)));


            Console.WriteLine($"\nsource_q_coords :\n{source_q_coords}");
            //var target_q_coords = np.concatenate((target_coords, np.zeros((1, n))));
            var target_q_coords = np.concatenate((target_coords, np.zeros(new shape(1, n), np.Int32)));

            Console.WriteLine($"\ntarget_q_coords :\n{target_q_coords}");

            //eq 40
            var b_scalar = _get_scalar(alpha_0, source_q_coords);
            Console.WriteLine($"\nb_scalar :\n{b_scalar}");

            //eq 41
            //var c_scalar = np.einsum("i->", temp);
            var c_scalar = np.sum(alpha_0, axis: -1);
            Console.WriteLine($"\nc_scalar :\n{c_scalar}");

            //eq 17
            var q0_w_matrix = _get_w_matrix(source_q_coords.T);
            Console.WriteLine($"\nq0_w_matrix :\n{q0_w_matrix}");

            //eq 16
            var qt_q_matrix = _get_q_matrix(target_q_coords.T);
            Console.WriteLine($"\nqt_q_matrix :\n{qt_q_matrix}");

            //eq 42
            var a_matrix = _get_abc_matrices(alpha_0, q0_w_matrix, qt_q_matrix);
            Console.WriteLine($"\na_matrix :\n{a_matrix}");

            //eq 43
            var b_matrix = _get_abc_matrices(alpha_0, qt_q_matrix);
            Console.WriteLine($"\nb_matrix :\n{b_matrix}");

            //eq 44
            var c_matrix = _get_abc_matrices(alpha_0, q0_w_matrix);
            Console.WriteLine($"\nc_matrix :\n{c_matrix}");


            double lambda_i = lambda_0;
            double i = 1;

            ndarray blc_matrix, d_matrix, beta_1, r_quat;
            (blc_matrix, d_matrix, beta_1, r_quat, lambda_i, i) = _get_solution(a_matrix,
                                                                              b_scalar,
                                                                              b_matrix,
                                                                              c_scalar,
                                                                              c_matrix,
                                                                              scale,
                                                                              lambda_i,
                                                                              i);

            Console.WriteLine($"\nblc_matrix :\n{blc_matrix}");
            Console.WriteLine($"\nd_matrix :\n{d_matrix}");
            Console.WriteLine($"\nbeta_1 : {beta_1}");
            Console.WriteLine($"\nr_quat :\n{r_quat}");
            Console.WriteLine($"\nlambda_i : {lambda_i}");
            Console.WriteLine($"\ni : {i}");

            //eq 25
            var r_matrix = _get_r_matrix(r_quat);
            Console.WriteLine($"\nr_matrix :\n{r_matrix}");

            //eq 58
            var s_quat = _get_s_quat(c_scalar, blc_matrix, r_quat);
            Console.WriteLine($"\ns_quat :\n{s_quat}");

            //eq 28, 26
            var t_vector = np.array(_get_t_vector(r_quat, s_quat)).reshape(3, 1);
            Console.WriteLine($"\nt_vector :\n{t_vector}");

            //Outputs
            var m = lambda_i;
            var r = r_matrix;
            var t = t_vector;

            //Verification
            var target_computed_points = (m * np.matmul(r, source_coords) + t).T;
            Console.WriteLine($"\ntarget_computed_points :\n{target_computed_points}");

        }




    }
}
