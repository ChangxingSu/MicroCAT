import argparse
import concurrent.futures
import os

# 定义过滤函数，用于处理每行数据
def filter_data(line):
    columns = line.strip().split("\t")
    if len(columns) > 0 and columns[0] == "C":
        return line
    return None

def extract(input_file, output_file, cores):
    # 使用多线程或多进程加快处理速度
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        # 获取系统可用核数
        num_cores = min(os.cpu_count(), cores)
        
        # 创建线程池或进程池
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:
            # 使用map方法将过滤函数应用于每一行数据
            results = executor.map(filter_data, infile)
            
            # 将满足条件的数据写入新表文件
            for result in results:
                if result is not None:
                    outfile.write(result)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter data based on the first column.")
    parser.add_argument("--input_file", help="Input TSV file name.")
    parser.add_argument("--output_file", help="Output TSV file name.")
    parser.add_argument("--cores", type=int, default=8,
                        help="Number of CPU cores to use for parallel processing.")
    args = parser.parse_args()
    print("run")
    extract(args.input_file, args.output_file, args.cores)
    print("fiish")

# 使用awk过滤数据并保存到新文件

