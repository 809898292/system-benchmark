import csv
import random
import time
from datetime import datetime, timedelta

# 配置参数
OUTPUT_FILE = "stock_data.csv"  # 输出文件名
TARGET_SIZE_MB = 100            # 目标文件大小 (MB)
STOCK_SYMBOLS = ["AAPL", "GOOGL", "MSFT", "AMZN", "TSLA", "NVDA", "META", "NFLX"]  # 股票代码列表

# 生成随机股票数据
def generate_stock_data(start_price, days=365):
    current_price = start_price
    data = []
    current_date = datetime.now() - timedelta(days=days)
    
    for _ in range(days):
        # 模拟价格波动 (-5% 到 +5%)
        change_percent = random.uniform(-0.05, 0.05)
        close_price = current_price * (1 + change_percent)
        high_price = max(current_price, close_price) * random.uniform(1.0, 1.03)
        low_price = min(current_price, close_price) * random.uniform(0.97, 1.0)
        open_price = current_price
        
        # 模拟成交量 (1M 到 50M 之间)
        volume = random.randint(1_000_000, 50_000_000)
        
        # 添加到数据
        data.append([
            current_date.strftime("%Y-%m-%d"),  # 日期
            random.choice(STOCK_SYMBOLS),      # 股票代码
            f"{open_price:.2f}",               # 开盘价
            f"{high_price:.2f}",               # 最高价
            f"{low_price:.2f}",               # 最低价
            f"{close_price:.2f}",              # 收盘价
            str(volume)                        # 成交量
        ])
        
        current_price = close_price
        current_date += timedelta(days=1)
    
    return data

# 写入CSV文件
def write_to_csv(filename, data):
    with open(filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Date", "Symbol", "Open", "High", "Low", "Close", "Volume"])  # 表头
        writer.writerows(data)

# 生成足够大的数据 (直到文件大小接近 100MB)
def generate_large_stock_data():
    data = []
    file_size = 0
    chunk_size = 100_000  # 每次生成 10 万行
    
    print("Generating stock data...")
    
    while file_size < TARGET_SIZE_MB * 1024 * 1024:
        chunk = []
        for _ in range(chunk_size):
            start_price = random.uniform(50, 5000)  # 随机初始价格 (50-5000)
            chunk.extend(generate_stock_data(start_price, days=1))  # 生成 1 天数据
        
        data.extend(chunk)
        
        # 检查当前文件大小
        with open(OUTPUT_FILE, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["Date", "Symbol", "Open", "High", "Low", "Close", "Volume"])
            writer.writerows(data)
        
        file_size = os.path.getsize(OUTPUT_FILE)
        print(f"Generated: {file_size / (1024 * 1024):.2f} MB", end="\r")
    
    print(f"\nDone! Output file: {OUTPUT_FILE} ({file_size / (1024 * 1024):.2f} MB)")

if __name__ == "__main__":
    import os
    generate_large_stock_data()
