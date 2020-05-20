using System;
using System.Collections.Generic;
using System.Linq;

public class FMM
{
    enum Label { Far, Considered, Accepted };

    public class PriorityQueue<T1, T2>
    where T1 : IComparable<T1>
    where T2 : IComparable<T2>
    {
        private SortedDictionary<T1, List<T2>> _queues_internal;

        public PriorityQueue()
        {
            _queues_internal = new SortedDictionary<T1, List<T2>>();
        }

        public int Count() { return _queues_internal.Count; }

        public void PushToQueue(T2 index, T1 dist)
        {
            List<T2> indices;
            if (_queues_internal.TryGetValue(dist, out indices))
            {
                indices.Add(index);
            }
            else
            {
                indices = new List<T2>();
                indices.Add(index);
                _queues_internal.Add(dist, indices);
            }

        }

        public void RemoveFromQueue(T2 index, T1 dist)
        {
            List<T2> indices;
            if (_queues_internal.TryGetValue(dist, out indices))
            {
                var location = indices.FindIndex(x => index.Equals(x));
                if (-1 != location)
                {
                    indices.RemoveAt(location);

                    if (0 == indices.Count)
                    {
                        if (!_queues_internal.Remove(dist))
                        {
                            throw new Exception("ERROR: cannot delete, why?");
                        }
                    }
                    return;
                }
            }

            throw new Exception("ERROR: cannot find this node from queue - node #" + index);
        }

        public T2 Dequeue()
        {
            var KV = _queues_internal.First();

            List<T2> indices = KV.Value;
            T2 index = indices[0];

            indices.RemoveAt(0);
            if (0 == indices.Count)
            {
                if (!_queues_internal.Remove(KV.Key))
                {
                    throw new Exception("ERROR: cannot delete, why?");
                }
            }

            return index;
        }
    }

    public static double[] Compute(bool[] inVoxels, int voX, int voY, int voZ, int[] voxel_indices_init, double[] distances_init)
    {
        // set initial status
        Func<int, int, int, int> GetIndex = (i, j, k) => { return i + j * voX + k * voX * voY; };

        ////////////////////////////////////////////////////////////
        // 1. Assign every node x_i the value of U_i = +inf and label them as FAR,
        //    for all nodes x_i \in \del\omega set U_i=0 and label x_i as ACCEPTED
        ////////////////////////////////////////////////////////////

        // initialize all nodes as far node {Far, +inf}
        double sdf_inf = double.MaxValue;
        int voXYZ = voX * voY * voZ;
        Label[] label = Enumerable.Repeat(Label.Far, voXYZ).ToArray();
        double[] sdf = Enumerable.Repeat(sdf_inf, voXYZ).ToArray();

        for (int idx = 0; idx < inVoxels.Length; idx++)
        {
            // do not recalculate distance for this node
            if (!inVoxels[idx])
            {
                label[idx] = Label.Accepted;
                sdf[idx] = float.MaxValue;
            }
        }

        // initialize priority queue
        PriorityQueue<double, int> pqueue = new PriorityQueue<double, int>();

        for (int i = 0; i < voxel_indices_init.Length; i++)
        {
            int v_index = voxel_indices_init[i];
            double v_distance = distances_init[i];

            label[v_index] = Label.Accepted; // still not sure
            sdf[v_index] = v_distance;

            pqueue.PushToQueue(v_index, v_distance);
        }

        ////////////////////////////////////////////////////////////
        // 2. for every far node x_i, use the "Eikonal update formula" to calculate a new value for U~.
        //    if U~ < U_i then set U_i = U~ and label x_i as CONSIDERED.
        ////////////////////////////////////////////////////////////
        Func<int, int, int, double> EikonalUpdate3D = (i, j, k) =>
        {
            double sdf_011 = (i == 0) ? sdf_inf : sdf[GetIndex(i - 1, j, k)];
            double sdf_211 = (i == voX - 1) ? sdf_inf : sdf[GetIndex(i + 1, j, k)];
            double sdf_101 = (j == 0) ? sdf_inf : sdf[GetIndex(i, j - 1, k)];
            double sdf_121 = (j == voY - 1) ? sdf_inf : sdf[GetIndex(i, j + 1, k)];
            double sdf_110 = (k == 0) ? sdf_inf : sdf[GetIndex(i, j, k - 1)];
            double sdf_112 = (k == voZ - 1) ? sdf_inf : sdf[GetIndex(i, j, k + 1)];

            List<double> phi = new List<double> { Math.Min(sdf_011, sdf_211), Math.Min(sdf_101, sdf_121), Math.Min(sdf_110, sdf_112) };
            phi.Sort();

            // from "Fluid Simulation for Computer Graphics," in 56p.
            //       written by Robert Bridson (ISBN 9781482232837)
            const double dx = 1.0;

            // Try just the closest neighbor
            double d = phi[0] + dx;
            if (d > phi[1])
            {
                // Try the two closest neighbors
                double sum = phi[0] + phi[1];
                d = 1.0 / 2.0 * (sum + Math.Sqrt(2.0f * dx * dx - (phi[1] - phi[0]) * (phi[1] - phi[0])));
                if (d > phi[2])
                {
                    // Use all three neighbors
                    sum = phi[0] + phi[1] + phi[2];
                    d = 1.0f / 3.0f * (sum + Math.Sqrt(Math.Max(0.0f, sum * sum - 3.0f * (phi[0] * phi[0] + phi[1] * phi[1] + phi[2] * phi[2] - dx * dx))));
                }
            }

            return d;
        };

        {
            double[] sdf_temp = (double[])sdf.Clone();
            Label[] label_temp = (Label[])label.Clone();

            for (int k = 0; k < voZ; k++)
                for (int j = 0; j < voY; j++)
                    for (int i = 0; i < voX; i++)
                    {
                        int idx = GetIndex(i, j, k);
                        if (Label.Far == label[idx])
                        {
                            double U_bar = EikonalUpdate3D(i, j, k);

                            if (U_bar < sdf[idx])
                            {
                                // set node as "CONSIDERED" node {U_bar, Considered}
                                if (Label.Considered == label[idx])
                                {
                                    pqueue.RemoveFromQueue(idx, sdf[idx]);
                                }

                                sdf_temp[idx] = U_bar;
                                label_temp[idx] = Label.Considered;
                                pqueue.PushToQueue(idx, U_bar);
                            }
                        }
                    }

            sdf = sdf_temp;
            label = label_temp;
        }

        ////////////////////////////////////////////////////////////
        // main loop
        ////////////////////////////////////////////////////////////
        while (0 < pqueue.Count())
        {
            ////////////////////////////////////////////////////////////
            // 3. let x~ be the CONSIDERED node with the smallest value U.
            //    label x~ as ACCEPTED.
            ////////////////////////////////////////////////////////////
            int id_u = pqueue.Dequeue();
            label[id_u] = Label.Accepted;

            int i1 = id_u % voX;
            int j1 = ((id_u - i1) / voX) % voY;
            int k1 = ((id_u - i1) / voX) / voY;
            double dist_u = sdf[id_u];

            ////////////////////////////////////////////////////////////
            // 4. for every neightbor x_i of x~ that is not-accepted,
            //    calculate a tentative value U~
            ////////////////////////////////////////////////////////////

            // 6 neighbors from 3D regular grid
            int[] di = { +1, -1, +0, +0, +0, +0 };
            int[] dj = { +0, +0, +1, -1, +0, +0 };
            int[] dk = { +0, +0, +0, +0, +1, -1 };

            for (int n = 0; n < 6; n++)
            {
                int i = i1 + di[n];
                int j = j1 + dj[n];
                int k = k1 + dk[n];

                // check validity of array index
                if (i < 0 || voX <= i || j < 0 || voY <= j || k < 0 || voZ <= k)
                    continue;

                ////////////////////////////////////////////////////////////
                // 5. if U~ < U_i then set U_i = U~.
                //    if x_i was labeled as FAR, update the label to CONSIDERED.
                ////////////////////////////////////////////////////////////
                int idx = GetIndex(i, j, k);

                if (Label.Accepted != label[idx])
                {
                    double U_bar = EikonalUpdate3D(i, j, k);
                    if (U_bar < sdf[idx])
                    {
                        // if it already exists in queue, find and remove it !
                        if (Label.Considered == label[idx])
                        {
                            pqueue.RemoveFromQueue(idx, sdf[idx]);
                            pqueue.PushToQueue(idx, U_bar);
                        }
                        else
                        {
                            pqueue.PushToQueue(idx, U_bar);
                        }

                        sdf[idx] = U_bar;
                        label[idx] = Label.Considered;
                    }
                }
            }
        }

        return sdf;
    }
}
