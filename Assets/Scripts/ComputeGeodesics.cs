using System.Collections.Generic;
using UnityEngine;

using System;
using System.Linq;
using System.Diagnostics;
using System.Runtime.CompilerServices;

public class ComputeGeodesics : MonoBehaviour
{
    #region PUBLIC_MEMBERS

    [Header("This option will have no effect during runtime.")]
    public string filepath = "Kitten_256x256.png";

    #endregion // PUBLIC_MEMBERS



    #region PROTECTED_MEMBERS

    Texture2D tex2D = null;

    int i0, j0;
    int W, H;
    bool[] area = null;

    #endregion // PROTECTED_MEMBERS



    #region MONO_BEHAVIOUR

    void Start()
    {
        LoadImage2D(filepath);
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

        {
            stopwatch.Start();
            double[] phi = FMM.Compute(area, W, H, 1, new int[1] { GetVoxelIndex(i0, j0) }, new double[1] { 0.0f } );
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

    void LoadImage2D(string filepath)
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
        for (int j1 = 0; j1 < H; j1++)
        {
            for (int i1 = 0; i1 < W; i1++)
            {
                int voxid_this = GetVoxelIndex(i1, j1);
                byte a = pixels32[voxid_this].a;

                if (a < 128)
                {
                    area[voxid_this] = true;
                }
            }
        }
    }

    #endregion // PRECOMPUTATION



    #region INLINE_FUNCTIONS

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private int GetVoxelIndex(int i, int j) => i + j * W;

    #endregion // INLINE_FUNCTIONS
}
