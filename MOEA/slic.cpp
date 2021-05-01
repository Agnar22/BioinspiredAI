#include "slic.h"
#include <iostream>
#include <math.h>
#include <queue>

void show_segmented_image(std::vector<int> closest_cluster, std::vector<slic::RGBCluster> &clusters, int w, int h) {
    cv::Mat img(h, w, CV_8UC3, cv::Scalar(0, 0, 0));
    for (int i=0; i<h*w; ++i) {
        int x = i%w;
        int y = i/w;
        slic::RGBCluster curr_clust = clusters[closest_cluster[i]];

        for (int clr_pos=0; clr_pos<3; ++clr_pos) {
            img.at<cv::Vec3b>(y, x)[clr_pos] = (unsigned char)(curr_clust.rgb[clr_pos]);
        }
        //std::cout << std::round(curr_clust.y) << " " <<  std::round(curr_clust.x) << std::endl;
        img.at<cv::Vec3b>(std::round(curr_clust.y), std::round(curr_clust.x)) = cv::Vec3b{0, 0, 255};
    }
    std::string name="Image";
    cv::namedWindow(name);
    cv::imshow(name, img);
    cv::waitKey(0);
}

namespace slic {
    std::vector<int> slic(cv::Mat &img, int num_segs) {
        std::vector<std::vector<cv::Point3_<uchar>>> img_c(img.rows, std::vector<cv::Point3_<uchar>>(img.cols));
        cv::Mat destination;
        cv::cvtColor(img, destination, cv::COLOR_RGB2Lab);
        int step = destination.step;
        int channels = destination.channels();
        for (int i = 0; i < destination.rows; i++) {
            for (int j = 0; j < destination.cols; j++) {
                cv::Point3_<uchar> pixelData;
                //L*: 0-255 (elsewhere is represented by 0 to 100)
                pixelData.x = (uchar)((double)(destination.data[step*i + channels*j + 0])/(255.0/100.0));
                //a*: 0-255 (elsewhere is represented by -127 to 127)
                pixelData.y = destination.data[step*i + channels*j + 1];
                //b*: 0-255 (elsewhere is represented by -127 to 127)
                pixelData.z = destination.data[step*i + channels*j + 2];
                std::cout << "bgr" << (int)img.at<cv::Vec3b>(i, j)[0] << " " << (int)img.at<cv::Vec3b>(i, j)[1] << " " << (int)img.at<cv::Vec3b>(i, j)[2] << std::endl;
                std::cout << "cielab" << (int)pixelData.x << " " << (int)pixelData.y << " " << (int)pixelData.z << std::endl;
                img_c[i][j] = pixelData;
            }
        }
        // Init clusters
        //std::cout << "start" << std::endl;
        int w = img.cols;
        int h = img.rows;

        double ratio = (double)(num_segs)/(double)(w*h);
        int num_cols = std::round(sqrt(ratio)*w);
        int num_rows = std::round(sqrt(ratio)*h);

        int s = std::round(((double)(w)/(num_cols+1)+(double)(h)/(num_rows+1))/2);
        int pad_left = (w - (num_cols - 1) * s)/2;
        int pad_top = (h - (num_rows - 1) * s)/2;
        std::cout << "s " << s << " pad_left " << pad_left << " pad_top " << pad_top << std::endl;
        std::cout << sqrt((w*h)/(num_rows*num_cols)) << std::endl;

        /*
        std::cout << w << " " << h << " " << ratio << std::endl;
        std::cout << num_cols << " " << num_rows << " " << num_segs << std::endl;
        std::cout << s << " " << pad_left << " " << pad_top << std::endl;
        */

        // Initialize clusters.
        //std::cout << "end" << std::endl;
        int num_clusters = num_cols*num_rows;
        std::vector<RGBCluster> clusters;
        for (int i=0; i<num_clusters; ++i) {
            //std::cout << i << std::endl;
            int clust_x = i%num_cols;
            int clust_y = i/num_cols;
            int y = clust_y==0 ? pad_top : pad_top+clust_y*s;
            int x = clust_x==0 ? pad_left : pad_left+clust_x*s;
            double r=0, g=0, b=0;
            for (int diff_x=-1; diff_x<=1; ++diff_x) {
                for (int diff_y=-1; diff_y<=1; ++diff_y) {
                    /*
                    r += (double)(img.at<cv::Vec3b>(y, x)[0]);
                    g += (double)(img.at<cv::Vec3b>(y, x)[1]);
                    b += (double)(img.at<cv::Vec3b>(y, x)[2]);
                    */
                    r += (double)img_c[y][x].x;
                    g += (double)img_c[y][x].y;
                    b += (double)img_c[y][x].z;

                }
            }
            //std::cout << "to_add" << std::endl;
            clusters.push_back(RGBCluster(r/9.0, g/9.0, b/9.0, x, y));

            //std::cout << r << " " << g << " " << b << " " << clusters[i].rgb[0] << " " << clusters[i].rgb[1] << " "<< clusters[i].rgb[2] << std::endl;
        }
        //std::cout << "end_loop" << std::endl;

        std::vector<int> pxl_cluster(w*h, -1);
        for (int iterations=0; iterations<100; ++iterations) {
            pxl_cluster.clear();
            pxl_cluster.resize(w*h, -1);
            std::vector<std::vector<double>> pxl_cluster_dist(w*h);
            std::vector<std::vector<double>> pxl_closest_clusters(w*h);
            for (int i=0; i<num_clusters; ++i) {
                /*
                int clust_x = i%num_cols;
                int clust_y = i/num_cols;
                int min_row = clust_y==0 ? 0 : pad_top+(clust_y-1)*s;
                int max_row = clust_y==num_rows-1 ? h : pad_top+(clust_y+1)*s;
                int min_col = clust_x==0 ? 0 : pad_left+(clust_x-1)*s;
                int max_col = clust_x==num_cols-1 ? w : pad_left+(clust_x+1)*s;
                */
                int min_row = std::max((int)(clusters[i].y-1*s), 0);
                int max_row = std::min((int)(clusters[i].y+1*s), h);
                int min_col = std::max((int)(clusters[i].x-1*s), 0);
                int max_col = std::min((int)(clusters[i].x+1*s), w);
                //std::cout << min_row << " " << max_row << " " << min_col << " " << max_col << std::endl;
                for (int row=min_row; row<max_row; ++row) {
                    for (int col=min_col; col<max_col; ++col) {
                        int pxl = row*w+col;
                        //double clr_dist = euc_dist(img.at<cv::Vec3b>(row, col), clusters[i].rgb);
                        double clr_dist = euc_dist(img_c[row][col], clusters[i].rgb);
                        double pos_dist = euc_dist((double)(col), (double)(row), clusters[i].x, clusters[i].y);
                        double dist = clr_dist + ((double)(10)/(double)(s)) * pos_dist;
                        pxl_cluster_dist[pxl].push_back(dist);
                        pxl_closest_clusters[pxl].push_back(i);
                    }
                }
            }
            for (int row=0; row<h; ++row) {
                for (int col=0; col<w; ++col) {
                    int pxl = row*w+col;
                    if (pxl_closest_clusters[pxl].size()>0)
                        continue;
                    std::cout << "not found" << std::endl;
                    for (int i=0; i<num_clusters; ++i) {
                        //double clr_dist = euc_dist(img.at<cv::Vec3b>(row, col), clusters[i].rgb);
                        double clr_dist = euc_dist(img_c[row][col], clusters[i].rgb);
                        double pos_dist = euc_dist((double)(col), (double)(row), clusters[i].x, clusters[i].y);
                        double dist = clr_dist + ((double)(10)/(double)(s)) * pos_dist;
                        pxl_cluster_dist[pxl].push_back(dist);
                        pxl_closest_clusters[pxl].push_back(i);
                    }
                }
            }
            std::vector<RGBCluster> next_clusters(num_clusters, RGBCluster());
            for (int row=0; row<h; ++row) {
                for (int col=0; col<w; ++col) {
                    int pxl = row*w+col;
                    int closest_cluster_pos = std::min_element(pxl_cluster_dist[pxl].begin(), pxl_cluster_dist[pxl].end())-pxl_cluster_dist[pxl].begin();
                    int closest_cluster = pxl_closest_clusters[pxl][closest_cluster_pos];


                    //std::cout << (int)img.at<cv::Vec3b>(row, col)[0] << " " << (int)img.at<cv::Vec3b>(row, col)[1] << " " << (int)img.at<cv::Vec3b>(row, col)[2] << std::endl;
                    next_clusters[closest_cluster].add_pxl(col, row, img_c[row][col]);//img.at<cv::Vec3b>(row, col));
                    pxl_cluster[pxl]=closest_cluster;
                    // Compute residual errors.
                }
            }
            show_segmented_image(pxl_cluster, next_clusters, w, h);
            clusters = next_clusters;
        }
       // TODO: Encforce connectivity.
       //return std::vector<int>(3, 0);
       //pxl_cluster = enforce_connectivity(pxl_cluster, clusters, w, h);
       show_segmented_image(pxl_cluster, clusters, w, h);
       return pxl_cluster;
    }

    std::vector<int> get_neighbours(int pos, int w, int h) {
        std::vector<int> neighbours;
        if (pos%w>0)
            neighbours.push_back(pos-1);
        if (pos%w<w-1)
            neighbours.push_back(pos+1);
        if (pos/w>0)
            neighbours.push_back(pos-w);
        if (pos/w<h-1)
            neighbours.push_back(pos+w);
        return neighbours;
    }

    std::vector<int> enforce_connectivity(std::vector<int> roots, std::vector<RGBCluster> clusters, int w, int h) {
        std::vector<bool> visited(roots.size(), false);
        std::queue<int> potential_loose;
        for (int cur_cluster=0; cur_cluster<clusters.size(); ++cur_cluster) {
            int x = std::round(clusters[cur_cluster].x);
            int y = std::round(clusters[cur_cluster].y);
            std::queue<int> q;
            q.push(y*w+x);
            visited[y*w+x] = true;
            while (!q.empty()) {
                int front = q.front();
                q.pop();
                for (int neighbour:get_neighbours(front, w, h)) {
                    if (roots[neighbour]==cur_cluster && !visited[neighbour]) {
                        visited[neighbour] = true;
                        q.push(neighbour);
                    } else if (roots[neighbour]!=cur_cluster && !visited[neighbour]) {
                        potential_loose.push(neighbour);
                    }

                }
            }
            //std::cout << "a" << std::endl;
        }

        // BFS that counts number of each neighbour.
        while (!potential_loose.empty()) {
            int start = potential_loose.front();
            potential_loose.pop();
            if (visited[start])
                continue;

            int cur_cluster = roots[start];
            std::queue<int> q;
            q.push(start);
            std::vector<int> others(clusters.size(), 0);
            while (!q.empty()) {
                //std::cout << "b" << std::endl;
                int front = q.front();
                q.pop();
                for (int neighbour:get_neighbours(front, w, h)) {
                    if (roots[neighbour]==cur_cluster && !visited[neighbour]) {
                        visited[neighbour] = true;
                        q.push(neighbour);
                    } else if (roots[neighbour]!=cur_cluster && visited[neighbour]) {
                        others[roots[neighbour]]++; 
                    }
                }
            }
            int most_frequent_neighbour = std::min_element(others.begin(), others.end())-others.begin();
            if (most_frequent_neighbour==roots[start])
                continue;
            std::queue<int> to_convert;
            to_convert.push(start);
            roots[start] = most_frequent_neighbour;
            //std::cout << "c" << std::endl;
            while (!to_convert.empty()) {
                int front = to_convert.front();
                to_convert.pop();
                for (int neighbour:get_neighbours(front, w, h)) {
                    if (roots[neighbour]==cur_cluster && visited[neighbour]) {
                        roots[neighbour] = most_frequent_neighbour;
                        to_convert.push(neighbour);
                    }
                }
            }
            //std::cout << "out" << std::endl;
        }
        return roots;
    }
}