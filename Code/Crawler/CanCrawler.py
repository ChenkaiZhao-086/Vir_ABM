import pandas as pd
import os
import re
from curl_cffi import requests as cffi_requests
from pyquery import PyQuery as pq
import time
import random
from datetime import datetime
from urllib.parse import urljoin


# 配置
BASE_URL = "https://www.canada.ca"
START_URL = "https://www.canada.ca/en/public-health/services/surveillance/respiratory-virus-detections-canada.html"
OUTPUT_DIR = "canada"


def get_html(url):
    """封装请求，增加简单的重试机制"""
    try:
        # impersonate="chrome110" 会让你的 TLS 指纹和 Chrome 110 完全一致
        resp = cffi_requests.get(url, impersonate="chrome110", timeout=10)

        if resp.status_code == 200:
            return resp.text
        else:
            print(f"Status Code: {resp.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching {url}: {e}")
        return None


def parse_date_from_url(url):
    """
    从URL中提取日期并转换为 YYYY-MM-DD 格式
    URL 示例: .../week-18-ending-may-7-2022.html
    """
    # 匹配 ending-month-day-year 模式
    match = re.search(r"ending-([a-zA-Z]+-\d{1,2}-\d{4})", url)
    if match:
        date_str = match.group(1)
        try:
            # %B 匹配全名月份 (e.g., May, November)
            dt = datetime.strptime(date_str, "%B-%d-%Y")
            return dt.strftime("%Y-%m-%d")
        except ValueError:
            pass
    return None


def get_all_report_links():
    """
    1. 访问主页
    2. 获取 Archives 下的所有年份页面
    3. 获取所有页面下的 Weekly Reports 链接
    """
    print("正在获取报告列表...")
    main_html = get_html(START_URL)
    if not main_html:
        return []

    doc = pq(main_html)

    # 1. 收集所有需要扫描的页面（主页 + Archives 中的年份页）
    pages_to_scan = [START_URL]

    # 在主页寻找 Archives 部分的链接 (根据提供的文本，链接在 "Archives" 标题下)
    # 假设 Archives 链接包含 'respiratory-virus-detections-canada' 且是年份格式
    archive_links = doc('a[href*="respiratory-virus-detections-canada"]').items()
    for link in archive_links:
        href = link.attr("href")
        # 简单的逻辑判断是否为年份归档页 (包含数字年份)
        if re.search(r"\d{4}-\d{4}", href):
            full_url = urljoin(BASE_URL, href)
            if full_url not in pages_to_scan:
                pages_to_scan.append(full_url)

    # 2. 扫描所有页面获取具体的周报链接
    weekly_report_urls = set()  # 使用集合去重

    for page_url in pages_to_scan:
        print(f"Scanning: {page_url}")
        page_html = get_html(page_url)
        if page_html:
            sub_doc = pq(page_html)
            # 提取所有包含 'week' 的链接
            links = sub_doc('a[href*="week"]').items()
            for link in links:
                href = link.attr("href")
                full_url = urljoin(BASE_URL, href)
                weekly_report_urls.add(full_url)

        time.sleep(random.randint(1, 2))

    return list(weekly_report_urls)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 获取所有周报链接
    all_urls = get_all_report_links()
    print(f"总共找到 {len(all_urls)} 个周报链接。")

    for url in all_urls:
        # 1. 提前解析日期，确定文件名
        date_filename = parse_date_from_url(url)

        # 如果无法解析日期，使用原始文件名作为备选
        if not date_filename:
            print(f"无法解析日期: {url}，跳过或使用默认名")
            continue

        file_path = os.path.join(OUTPUT_DIR, f"{date_filename}.csv")

        # 2. 断点续传：如果文件已存在，跳过
        if os.path.exists(file_path):
            print(f"文件已存在，跳过: {file_path}")
            continue

        print(f"正在下载: {date_filename} ...")

        # 3. 下载并解析表格
        page_html = get_html(url)
        if not page_html:
            continue

        try:
            # 使用 pandas 直接寻找包含 "Respiratory virus detections" 的表格
            # match 参数可以确保我们抓取的是正确的表，而不是页面的布局表格
            dfs = pd.read_html(page_html, match="Respiratory virus detections")

            if dfs:
                target_df = dfs[0]
                target_df.to_csv(file_path, index=False)
                print(f"成功保存: {file_path}")
            else:
                print(f"未找到表格: {url}")

        except Exception as e:
            print(f"解析表格失败 {url}: {e}")

        # 礼貌爬取
        time.sleep(random.randint(1, 3))


if __name__ == "__main__":
    main()
