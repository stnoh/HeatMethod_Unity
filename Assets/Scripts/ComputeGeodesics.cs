using System.Collections.Generic;
using UnityEngine;

using System;
using System.Linq;
using System.Diagnostics;
using System.Runtime.CompilerServices;

// sparse matrix module (LGPL)
using CSparse;
using CSparse.Double;
using CSparse.Double.Factorization;
using CSparse.Storage;

public class ComputeGeodesics : MonoBehaviour
{
    #region PUBLIC_MEMBERS

    public enum MethodType { FMM, HEAT };
    public MethodType method = MethodType.HEAT;

    [Header("This option will have no effect during runtime.")]
    public string filepath = "Kitten_256x256.png";
    public double dt = 32.0; // small time for heat transfer

    #endregion // PUBLIC_MEMBERS



    #region PROTECTED_MEMBERS

    Texture2D tex2D = null;

    int i0, j0;
    int W, H;
    bool[] area = null;

    // map between voxel and matrix index
    int area_cnt;
    Dictionary<int, int> Voxel2Matrix;
    Dictionary<int, int> Matrix2Voxel;
    Dictionary<int, List<int>> MatrixNeighbors;

    #endregion // PROTECTED_MEMBERS



    #region PRIVATE_MEMBERS

    private SparseCholesky solver_heat;
    private SparseCholesky solver_dist;

    #endregion // PRIVATE_MEMBERS



    #region MONO_BEHAVIOUR

    void Start()
    {
        var stopwatch = new Stopwatch();

        PrecomputeValid2D(filepath);

        // step I-1. & I-2.
        stopwatch.Start();
        PrecomputeHeatMatrix2D(dt);
        stopwatch.Stop();
        print(string.Format("precompute for heat equation: {0} [msec]", stopwatch.ElapsedMilliseconds));

        // step III-1. & III-2.
        stopwatch.Reset();
        stopwatch.Start();
        PrecomputeDistanceMatrix2D(W, H);
        stopwatch.Stop();
        print(string.Format("precompute for distance equation: {0} [msec]", stopwatch.ElapsedMilliseconds));
    }

    void Update()
    {
        // left clicked position
        if (Input.GetMouseButton(0))
        {
            float mousePosX = (Input.mousePosition.x) / Screen.width;
            float mousePosY = (Input.mousePosition.y) / Screen.height;
            i0 = (int)(W * mousePosX);
            j0 = (int)(H * mousePosY);

            double[] phi_dense = CalcDistance();
            VisualizeDistanceDense(phi_dense);
        }
    }

    #endregion // MONO_BEHAVIOUR



    #region RUNTIME

    double[] CalcDistance()
    {
        var stopwatch = new Stopwatch();

        if (MethodType.HEAT == method)
        {
            // step I-3. solve heat equation
            stopwatch.Start();
            double[] u_heat_sparse = SolveHeat();

            // step II. normalized & negated gradients of heat
            double[][] X_sparse = ComputeGrad2D(u_heat_sparse);

            // step III-3. solve distance
            double[] phi = SolveDistance(X_sparse);
            stopwatch.Stop();

            print(string.Format("solve distance: {0} [msec]", stopwatch.ElapsedMilliseconds));
            return phi;
        }
        else
        {
            stopwatch.Start();
            double[] phi = FMM.Compute(area, W, H, 1, new int[1] { GetVoxelIndex(i0, j0) }, new double[1] { 0.0f });
            stopwatch.Stop();

            // handling out-of-range
            for (int i = 0; i < W * H; i++)
            {
                if (!area[i]) phi[i] = 0.0;
            }

            print(string.Format("fast marching method: {0} [msec]", stopwatch.ElapsedMilliseconds));
            return phi;
        }
    }

    // STEP I-3.
    double[] SolveHeat()
    {
        double init_heat = (double)System.Math.Sqrt(W * H) / 64.0;

        double[] b_heat_sparse = new double[area_cnt];
        double[] u_heat_sparse = new double[area_cnt];

        int voxid_init = GetVoxelIndex(i0, j0);
        int matid_init;
        if (Voxel2Matrix.TryGetValue(voxid_init, out matid_init))
        {
            b_heat_sparse[matid_init] = init_heat; // the point of initial seed
        }

        solver_heat.Solve(b_heat_sparse, u_heat_sparse);

        return u_heat_sparse;
    }

    // STEP II.
    double[][] ComputeGrad2D(double[] u_heat_sparse)
    {
        double[][] X_sparse = new double[2][];
        X_sparse[0] = new double[area_cnt];
        X_sparse[1] = new double[area_cnt];

        for (int matid_this = 0; matid_this < area_cnt; matid_this++)
        {
            List<int> matid_neighbors;

            if (!MatrixNeighbors.TryGetValue(matid_this, out matid_neighbors)) continue;

            double u_heat_c = u_heat_sparse[matid_this];

            double u_heat_l = (matid_neighbors[0] == -1) ? u_heat_c : u_heat_sparse[matid_neighbors[0]];
            double u_heat_r = (matid_neighbors[1] == -1) ? u_heat_c : u_heat_sparse[matid_neighbors[1]];
            double u_heat_b = (matid_neighbors[2] == -1) ? u_heat_c : u_heat_sparse[matid_neighbors[2]];
            double u_heat_t = (matid_neighbors[3] == -1) ? u_heat_c : u_heat_sparse[matid_neighbors[3]];

            double dx = (double)(u_heat_r - u_heat_l);
            double dy = (double)(u_heat_t - u_heat_b);
            double mag = System.Math.Sqrt(dx * dx + dy * dy);

            X_sparse[0][matid_this] = -dx / mag;
            X_sparse[1][matid_this] = -dy / mag;
        }

        return X_sparse;
    }

    // STEP III-3.
    double[] SolveDistance(double[][] X_sparse)
    {
        ////////////////////////////////////////////////////////////
        // divergence of normalized gradient X
        // [CAUTION] NEGATED value to fit Cholesky decomposition
        ////////////////////////////////////////////////////////////
        double[] b_dist_sparse = new double[area_cnt];

        // change for loop based on sparse
        for (int matid_this = 0; matid_this < area_cnt; matid_this++)
        {
            List<int> matid_neighbor;
            if (!MatrixNeighbors.TryGetValue(matid_this, out matid_neighbor)) continue;

            double div_x = 0.0;
            double div_y = 0.0;

            bool bound_x0 = false;
            bool bound_x2 = false;
            if (matid_neighbor[0] != -1) div_x -= X_sparse[0][matid_neighbor[0]]; else bound_x0 = true;
            if (matid_neighbor[1] != -1) div_x += X_sparse[0][matid_neighbor[1]]; else bound_x2 = true;
            if (bound_x0 || bound_x2) div_x *= 2.0;

            bool bound_y0 = false;
            bool bound_y2 = false;
            if (matid_neighbor[2] != -1) div_y -= X_sparse[1][matid_neighbor[2]]; else bound_y0 = true;
            if (matid_neighbor[3] != -1) div_y += X_sparse[1][matid_neighbor[3]]; else bound_y2 = true;
            if (bound_y0 || bound_y2) div_y *= 2.0;

            b_dist_sparse[matid_this] = -0.5 * (div_x + div_y); // NEGATED!
        }

        ////////////////////////////////////////////////////////////
        // solve distance field
        ////////////////////////////////////////////////////////////
        double[] phi_sparse = new double[area_cnt];
        solver_dist.Solve(b_dist_sparse, phi_sparse);

        double phi_max = phi_sparse.Max();
        double phi_min = phi_sparse.Min();

        double[] phi_dense = new double[W * H];

        for (int matid_this = 0; matid_this < area_cnt; matid_this++)
        {
            int voxid_this;
            Matrix2Voxel.TryGetValue(matid_this, out voxid_this);

            phi_dense[voxid_this] = phi_sparse[matid_this] - phi_min;
        }
        return phi_dense;
    }

    void VisualizeDistanceDense(double[] phi_dense)
    {
        Color32[] pixels32 = tex2D.GetPixels32();

        double phi_max = phi_dense.Max();
        double phi_min = phi_dense.Min();

        for (int j1 = 0; j1 < H; j1++)
        {
            for (int i1 = 0; i1 < W; i1++)
            {
                int voxid_this = GetVoxelIndex(i1, j1);

                if (area[voxid_this])
                {
                    float val01 = (float)((phi_dense[voxid_this] - phi_min) / (phi_max - phi_min));

                    // "hot" colormap in MATLAB/matplotlib
                    float r = Mathf.Clamp(2.62518523f * (255.0f * val01 -   0.0f) + 10.60800f, 0.0f, 255.0f); // at   0/255
                    float g = Mathf.Clamp(2.62499573f * (255.0f * val01 -  94.0f) +  2.37524f, 0.0f, 255.0f); // at  94/255
                    float b = Mathf.Clamp(3.93750394f * (255.0f * val01 - 191.0f) +  2.99975f, 0.0f, 255.0f); // at 191/255

                    pixels32[voxid_this].r = (byte)r;
                    pixels32[voxid_this].g = (byte)g;
                    pixels32[voxid_this].b = (byte)b;
                }
            }
        }

        tex2D.SetPixels32(pixels32);
        tex2D.Apply();
    }

    #endregion // RUNTIME



    #region PRECOMPUTATION

    void PrecomputeValid2D(string filepath)
    {
        // load & apply texture
        byte[] data = System.IO.File.ReadAllBytes(filepath);

        tex2D = new Texture2D(1, 1);
        tex2D.LoadImage(data);
        Color32[] pixels32 = tex2D.GetPixels32();

        gameObject.GetComponent<Renderer>().material.mainTexture = tex2D;

        // create boolean array for quick checking
        W = tex2D.width;
        H = tex2D.height;
        area = Enumerable.Repeat<bool>(false, W * H).ToArray();

        // initialize array
        area_cnt = 0;
        Voxel2Matrix = new Dictionary<int, int>();
        Matrix2Voxel = new Dictionary<int, int>();
        for (int j1 = 0; j1 < H; j1++)
        {
            for (int i1 = 0; i1 < W; i1++)
            {
                int voxid_this = GetVoxelIndex(i1, j1);
                byte a = pixels32[voxid_this].a;

                if (a < 128)
                {
                    area[voxid_this] = true;

                    // map between voxel/matrix index
                    Voxel2Matrix.Add(voxid_this, area_cnt);
                    Matrix2Voxel.Add(area_cnt, voxid_this);
                    area_cnt += 1;
                }
            }
        }

        // neighbor information in sparse matrix
        MatrixNeighbors = new Dictionary<int, List<int>>();
        for (int j1 = 0; j1 < H; j1++)
        {
            for (int i1 = 0; i1 < W; i1++)
            {
                int voxid_this = GetVoxelIndex(i1, j1);
                int matid_this;

                if (Voxel2Matrix.TryGetValue(voxid_this, out matid_this))
                {
                    List<int> matid_neighbors = new List<int>();

                    Action<int, int> CheckNeighbor = (i, j) =>
                    {
                        int voxid_next = GetVoxelIndex(i, j);
                        int matid_next;
                        if (Voxel2Matrix.TryGetValue(voxid_next, out matid_next))
                        {
                            matid_neighbors.Add(matid_next);
                        }
                        else
                        {
                            matid_neighbors.Add(-1);
                        }
                    };

                    CheckNeighbor(i1 - 1, j1);
                    CheckNeighbor(i1 + 1, j1);
                    CheckNeighbor(i1, j1 - 1);
                    CheckNeighbor(i1, j1 + 1);

                    MatrixNeighbors.Add(matid_this, matid_neighbors);
                }
            }
        }
    }

    void PrecomputeHeatMatrix2D(double dt)
    {
        int N = area_cnt;

        ////////////////////////////////////////////////////////////
        // step I-1. build matrix A for heat diffusion
        ////////////////////////////////////////////////////////////
        var C_heat = new CoordinateStorage<double>(N, N, 5 * N); // banded matrix without exception

        for (int matid_this = 0; matid_this < N; matid_this++)
        {
            List<int> matid_neighbor;
            if (!MatrixNeighbors.TryGetValue(matid_this, out matid_neighbor)) continue;

            if (matid_neighbor[0] != -1) C_heat.At(matid_this, matid_neighbor[0], -1.0 * dt);
            if (matid_neighbor[1] != -1) C_heat.At(matid_this, matid_neighbor[1], -1.0 * dt);
            if (matid_neighbor[2] != -1) C_heat.At(matid_this, matid_neighbor[2], -1.0 * dt);
            if (matid_neighbor[3] != -1) C_heat.At(matid_this, matid_neighbor[3], -1.0 * dt);

            C_heat.At(matid_this, matid_this, 1.0 + 4.0 * dt);
        }

        var A_heat = Converter.ToCompressedColumnStorage(C_heat) as SparseMatrix;

        ////////////////////////////////////////////////////////////
        // step I-2. build matrix A for heat diffusion
        ////////////////////////////////////////////////////////////
        solver_heat = SparseCholesky.Create(A_heat, ColumnOrdering.MinimumDegreeAtPlusA);
    }

    void PrecomputeDistanceMatrix2D(int W, int H, bool show = true)
    {
        int N = area_cnt;

        ////////////////////////////////////////////////////////////
        // step III-1. build matrix A for distance computation
        // [CAUTION] this is NEGATED matrix to use Cholesky decomp.
        ////////////////////////////////////////////////////////////
        var C_dist = new CoordinateStorage<double>(N, N, 5 * N); // banded matrix without exception

        for (int matid_this = 0; matid_this < N; matid_this++)
        {
            List<int> matid_neighbor;
            if (!MatrixNeighbors.TryGetValue(matid_this, out matid_neighbor)) continue;

            double A_diag = -1e-6; // small value to use Cholesky decomposition

            Action<int> CheckNeighbor = (id) =>
            {
                int matid_next = matid_neighbor[id];
                if (matid_next != -1)
                {
                    A_diag += -1.0;
                    C_dist.At(matid_this, matid_next, -1.0);
                }
            };

            CheckNeighbor(0);
            CheckNeighbor(1);
            CheckNeighbor(2);
            CheckNeighbor(3);

            C_dist.At(matid_this, matid_this, -A_diag);
        }

        var A_dist = Converter.ToCompressedColumnStorage(C_dist) as SparseMatrix;

        ////////////////////////////////////////////////////////////
        // step III-2. build matrix A for heat diffusion
        ////////////////////////////////////////////////////////////
        solver_dist = SparseCholesky.Create(A_dist, ColumnOrdering.MinimumDegreeAtPlusA);
    }

    #endregion // PRECOMPUTATION



    #region INLINE_FUNCTIONS

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private int GetVoxelIndex(int i, int j) => i + j * W;

    #endregion // INLINE_FUNCTIONS
}
