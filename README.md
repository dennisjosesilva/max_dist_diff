# Compiling

1. Create a "build" directory and goes to this directory:
   
   ```shell
   $ mkdir build && cd build
   ```

2. Generate files for building:
   
   ```shell
   $ cmake .. 
   ```

3. build the programs:
   
   ```shell
   $ cmake --build .
   ```

4. run the programs:
   
   ```shell
   $ ./dtmt <image>
   $ ./dtmt_diff <image>
   $ ./max_dist_attr <in_image> <out_image>
   ```
