# Time: 2023-02-02
# Author: Haosen Luo
# Filename: ncbi_gene_table.py
# 注意事项, 需要保证环境中含有requests、faker、pandas等库
import os
import re
import time
import pandas as pd
import numpy as np
import requests
from faker import Faker
from configparser import ConfigParser
from requests.packages.urllib3.exceptions import InsecureRequestWarning

# from urllib3.util.retry import Retry
# from requests.adapters import HTTPAdapter

# 配置文件读取
BIN = os.path.dirname(__file__) + '/'
config = ConfigParser()
config.read(BIN + 'config.ini')

# 表格信息匹配
ncbi_genome_string = ['GRCh38.p14', 'T2T-CHM13v2.0', 'GRCh37.p13']
review_genome_string = ['hg38', 't2t', 'hg19']
review_genome_dict = dict(zip(ncbi_genome_string, review_genome_string))

# Unknown genome version
unknown_genome_version_extract_pattern = re.compile(r'.*\.(\w+)\s\(.*\)', re.I)

# Location提取表达式
# no complement
ncbi_genome_split_match_pattern = re.compile(r' \(', re.I)
ncbi_location_pattern = re.compile(r' ', re.I)
ncbi_location_no_complement_pos_pattern = re.compile(r'\((\d+)\.\.(\d+)\)')
# complement
ncbi_location_complement_pattern = re.compile(r'\(\d+\.\.\d+,.*?\)', re.I)
ncbi_location_complement_split_pos_pattern = re.compile(r'\((\d+)\.\.(\d+),.*?\)')

# 重连设置
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)


# retries = Retry(total=5,
#                 backoff_factor=0.1,
#                 status_forcelist=[500, 502, 503, 504])


class NcbiGeneParser(object):
    def __init__(self, ncbi_id=None, find_gene_name=None, save_html_dir=None):
        self.gene_id = str(ncbi_id) if isinstance(ncbi_id, int) else ncbi_id
        self.gene_name = find_gene_name
        if self.gene_id == '':
            raise Exception(f"gene_id is empty {self.gene_id}, please Re-enter the NCBI ID")

        gene_name_id = '-'.join([self.gene_name, self.gene_id])

        if save_html_dir:
            self.html_dir = save_html_dir if save_html_dir.endswith('/') else save_html_dir + '/'
            os.system(f"mkdir -p {self.html_dir}") if not os.path.exists(self.html_dir) else ...
            self.html_dir_gene_html_save_path = ''.join([self.html_dir, gene_name_id, '.html'])
        else:
            ...

        # url construction
        self.__gene_url = f"https://www.ncbi.nlm.nih.gov/gene/{self.gene_id}"

        self.hg19_transcript = None
        self.hg19_chr = None
        self.hg19_begin = None
        self.hg19_end = None
        self.hg38_transcript = None
        self.hg38_chr = None
        self.hg38_begin = None
        self.hg38_end = None
        self.t2t_transcript = None
        self.t2t_chr = None
        self.t2t_begin = None
        self.t2t_end = None
        # 额外未知
        self.unknown_transcript = None
        self.unknown_chr = None
        self.unknown_begin = None
        self.unknown_end = None

    @staticmethod
    def __read_html_content(url, max_retry=5, sleep_stop=30):
        # 调用
        faker = Faker()
        # 代理设置
        faker_agent_dict = {'user-agent': faker.user_agent()}
        # 重连设置
        switch = True
        now_retry = 0
        response_content = ''

        while switch:
            try:
                response = requests.get(url, headers=faker_agent_dict, timeout=100, stream=True, verify=False)
            except Exception:
                if now_retry == max_retry:
                    switch = False
                else:
                    time.sleep(sleep_stop)
                    now_retry += 1
            else:
                if response.status_code == 200:
                    switch = False
                    response_content = response.content
                    response.close()
                else:
                    time.sleep(sleep_stop)
                    now_retry += 1
        return response_content

    @property
    def save_html_content(self):
        url = self.__gene_url
        content = self.__read_html_content(url=url)
        html_save_path = self.html_dir_gene_html_save_path
        if html_save_path:
            with open(html_save_path, 'wb') as f:
                f.write(content)
        return None

    @staticmethod
    def __init_chromosome(chromosome):
        # 1 | '1' | 22 | '22' | xy | XY |
        if not re.search(re.compile(r'chr(\w+)', re.I), string=str(chromosome)):
            if isinstance(chromosome, int):
                return ''.join(['chr', str(chromosome)])
            elif isinstance(chromosome, str):
                if chromosome in [str(i) for i in np.arange(23)]:
                    return f'chr{chromosome}'
                elif chromosome in ['X', 'Y', 'x', 'y']:
                    return f"chr{chromosome.upper()}"
                else:
                    return ''
        else:
            return ''.join(['chr', str(re.split(re.compile(r'chr', re.I), string=str(chromosome).upper())[-1])])

    @staticmethod
    def __get_ref_genome_version_dict(assembly):
        ncbi_genome_version = re.search(unknown_genome_version_extract_pattern, assembly).groups()[0]
        if review_genome_dict.get(ncbi_genome_version, 'none') not in ['hg19', 'hg38', 't2t']:
            review_genome_dict[assembly] = ncbi_genome_version
        return ncbi_genome_version, review_genome_dict

    @staticmethod
    def __get_location_list(location_value):
        if location_value != '':
            # 分割单元格数据为列表
            info_list = re.split(ncbi_location_pattern, string=location_value)
            # 转录本
            transcript = info_list[0]
            # 针对location中是否有补充分情况汇总
            if 'complement' in location_value:
                complement_list = re.findall(ncbi_location_complement_pattern, string=location_value)
                if len(complement_list) == 1:
                    begin_pos, end_pos = re.search(ncbi_location_complement_split_pos_pattern,
                                                   string=complement_list[0]).groups()
                    return [(transcript, begin_pos, end_pos)]
                else:
                    transcript_begin_end_pos_list = []
                    transcript_begin_end_pos_list_append = transcript_begin_end_pos_list.append
                    for i in range(len(complement_list)):
                        begin_pos, end_pos = re.search(ncbi_location_complement_split_pos_pattern,
                                                       string=complement_list[i]).groups()
                        transcript_begin_end_pos_list_append((transcript, begin_pos, end_pos))
                    return transcript_begin_end_pos_list
            else:
                try:
                    begin_pos = re.search(ncbi_location_no_complement_pos_pattern, string=info_list[-1]).groups()[0]
                    end_pos = re.search(ncbi_location_no_complement_pos_pattern, string=info_list[-1]).groups()[-1]
                except Exception:
                    begin_pos, end_pos = '', ''
                return [(transcript, begin_pos, end_pos)]
        else:
            return [('', '', '')]

    @property
    def dataframe(self):
        df = pd.DataFrame()
        hg38_chr = ''
        hg38_transcript = ''
        hg38_begin = ''
        hg38_end = ''
        t2t_chr = ''
        t2t_transcript = ''
        t2t_begin = ''
        t2t_end = ''
        hg19_chr = ''
        hg19_transcript = ''
        hg19_begin = ''
        hg19_end = ''
        #
        ncbi_id = self.gene_id
        gene_name = self.gene_name
        url = self.__gene_url

        response_content = self.__read_html_content(url=url)
        if response_content:
            try:
                html_df = pd.read_html(response_content)[0]
            except Exception:
                html_df = pd.DataFrame(columns=['gene_name', 'ncbi_id',
                                                'hg38_chrom', 'hg38_transcript', 'hg38_begin', 'hg38_end',
                                                't2t_chrom', 't2t_transcript', 't2t_begin', 't2t_end',
                                                'hg19_chrom', 'hg19_transcript', 'hg19_begin', 'hg19_end', 'ncbi_url'],
                                       index=[0])
                html_df['gene_name'] = gene_name
                html_df['ncbi_id'] = ncbi_id
                html_df['ncbi_url'] = url
                return html_df
            else:
                if all(html_df.columns == ['Annotation release', 'Status', 'Assembly', 'Chr', 'Location']):
                    for idx, row in html_df.iterrows():
                        ncbi_genome_version = re.split(ncbi_genome_split_match_pattern, row['Assembly'])[0]
                        if review_genome_dict.get(ncbi_genome_version, 'none') == 'hg38':
                            if len(self.__get_location_list(row['Location'])) == 1:
                                self.hg38_transcript, self.hg38_begin, self.hg38_end =\
                                    self.__get_location_list(row['Location'])[0]
                                self.hg38_chr = self.__init_chromosome(row['Chr'])
                                hg38_chr, df.loc[0, 'hg38_chrom'] = self.hg38_chr, self.hg38_chr
                                hg38_transcript, df.loc[
                                    0, 'hg38_transcript'] = self.hg38_transcript, self.hg38_transcript
                                hg38_begin, df.loc[0, 'hg38_begin'] = self.hg38_begin, self.hg38_begin
                                hg38_end, df.loc[0, 'hg38_end'] = self.hg38_end, self.hg38_end
                            else:
                                for i in range(len(self.__get_location_list(row['Location']))):
                                    df.loc[i, 'hg38_chrom'] = self.__init_chromosome(row['Chr'])
                                    df.loc[i, 'hg38_transcript'] = self.__get_location_list(row['Location'])[i][0]
                                    df.loc[i, 'hg38_begin'] = self.__get_location_list(row['Location'])[i][1]
                                    df.loc[i, 'hg38_end'] = self.__get_location_list(row['Location'])[i][-1]
                        elif review_genome_dict.get(ncbi_genome_version, 'none') == 't2t':
                            if len(self.__get_location_list(row['Location'])) == 1:
                                self.t2t_transcript, self.t2t_begin, self.t2t_end =\
                                    self.__get_location_list(row['Location'])[0]
                                self.t2t_chr = self.__init_chromosome(row['Chr'])
                                t2t_chr, df.loc[0, 't2t_chrom'] = self.t2t_chr, self.t2t_chr
                                t2t_transcript, df.loc[
                                    0, 't2t_transcript'] = self.t2t_transcript, self.t2t_transcript
                                t2t_begin, df.loc[0, 't2t_begin'] = self.t2t_begin, self.t2t_begin
                                t2t_end, df.loc[0, 't2t_end'] = self.t2t_end, self.t2t_end
                            else:
                                for i in range(len(self.__get_location_list(row['Location']))):
                                    df.loc[i, 't2t_chrom'] = self.__init_chromosome(row['Chr'])
                                    df.loc[i, 't2t_transcript'] = self.__get_location_list(row['Location'])[i][0]
                                    df.loc[i, 't2t_begin'] = self.__get_location_list(row['Location'])[i][1]
                                    df.loc[i, 't2t_end'] = self.__get_location_list(row['Location'])[i][-1]
                        elif review_genome_dict.get(ncbi_genome_version, 'none') == 'hg19':
                            if len(self.__get_location_list(row['Location'])) == 1:
                                self.hg19_transcript, self.hg19_begin, self.hg19_end =\
                                    self.__get_location_list(row['Location'])[0]
                                self.hg19_chr = self.__init_chromosome(row['Chr'])
                                hg19_chr, df.loc[0, 'hg19_chrom'] = self.hg19_chr, self.hg19_chr
                                hg19_transcript, df.loc[
                                    0, 'hg19_transcript'] = self.hg19_transcript, self.hg19_transcript
                                hg19_begin, df.loc[0, 'hg19_begin'] = self.hg19_begin, self.hg19_begin
                                hg19_end, df.loc[0, 'hg19_end'] = self.hg19_end, self.hg19_end
                            else:
                                for i in range(len(self.__get_location_list(row['Location']))):
                                    df.loc[i, 'hg19_chrom'] = self.__init_chromosome(row['Chr'])
                                    df.loc[i, 'hg19_transcript'] = self.__get_location_list(row['Location'])[i][0]
                                    df.loc[i, 'hg19_begin'] = self.__get_location_list(row['Location'])[i][1]
                                    df.loc[i, 'hg19_end'] = self.__get_location_list(row['Location'])[i][-1]
                        else:
                            ncbi_genome_version, gene_version_dict = self.__get_ref_genome_version_dict(
                                row['Assembly'])
                            # ncbi_genome_version = re.search(unknown_genome_version_extract_pattern, row['Assembly']).groups()[0]
                            # if review_genome_dict.get(row['Assembly'], 'none') != ncbi_genome_version:
                            #     review_genome_dict[row['Assembly']] = ncbi_genome_version
                            self.unknown_chr = self.__init_chromosome(row['Chr'])
                            unknown_version_chr_str = '_'.join([ncbi_genome_version, 'chrom'])
                            unknown_version_transcript_str = '_'.join([ncbi_genome_version, 'transcript'])
                            unknown_version_begin_str = '_'.join([ncbi_genome_version, 'begin'])
                            unknown_version_end_str = '_'.join([ncbi_genome_version, 'end'])
                            #
                            df.loc[0, unknown_version_chr_str] = self.unknown_chr
                            if len(self.__get_location_list(row['Location'])) == 1:
                                self.unknown_transcript, self.unknown_begin, self.unknown_end =\
                                    self.__get_location_list(row['Location'])[0]
                                df.loc[0, unknown_version_transcript_str] = self.unknown_transcript
                                df.loc[0, unknown_version_begin_str] = self.unknown_begin
                                df.loc[0, unknown_version_end_str] = self.unknown_end
                            else:
                                for i in range(len(self.__get_location_list(row['Location']))):
                                    df.loc[i, unknown_version_transcript_str] =\
                                        self.__get_location_list(row['Location'])[i][0]
                                    df.loc[i, unknown_version_begin_str] = self.__get_location_list(row['Location'])[i][
                                        1]
                                    df.loc[i, unknown_version_end_str] = self.__get_location_list(row['Location'])[i][
                                        -1]
                else:
                    print(f'{self.gene_name}\t{self.gene_id}\t{self.__gene_url}\n')
                    df = pd.DataFrame(columns=['hg38_chrom', 'hg38_transcript', 'hg38_begin', 'hg38_end',
                                               't2t_chrom', 't2t_transcript', 't2t_begin', 't2t_end',
                                               'hg19_chrom', 'hg19_transcript', 'hg19_begin', 'hg19_end'],
                                      index=[0])

                df.fillna('-', inplace=True)
                for idx, line in df.iterrows():
                    # hg38 chrom
                    if 'hg38_chrom' in df.columns:
                        if line['hg38_chrom'] == '-':
                            df.loc[idx, 'hg38_chrom'] = hg38_chr
                    else:
                        df['hg38_chrom'] = np.nan
                    # hg38 transcript
                    if 'hg38_transcript' in df.columns:
                        if line['hg38_transcript'] == '-':
                            df.loc[idx, 'hg38_transcript'] = hg38_transcript
                    else:
                        df['hg38_transcript'] = np.nan
                    # hg38 begin
                    if 'hg38_begin' in df.columns:
                        if line['hg38_begin'] == '-':
                            df.loc[idx, 'hg38_begin'] = hg38_begin
                    else:
                        df['hg38_begin'] = np.nan
                    # hg38 end
                    if 'hg38_end' in df.columns:
                        if line['hg38_end'] == '-':
                            df.loc[idx, 'hg38_end'] = hg38_end
                    else:
                        df['hg38_end'] = np.nan
                    # t2t chrom
                    if 't2t_chrom' in df.columns:
                        if line['t2t_chrom'] == '-':
                            df.loc[idx, 't2t_chrom'] = t2t_chr
                    else:
                        df['t2t_chrom'] = np.nan
                    # t2t transcript
                    if 't2t_transcript' in df.columns:
                        if line['t2t_transcript'] == '-':
                            df.loc[idx, 't2t_transcript'] = t2t_transcript
                    else:
                        df['t2t_transcript'] = np.nan
                    # t2t begin
                    if 't2t_begin' in df.columns:
                        if line['t2t_begin'] == '-':
                            df.loc[idx, 't2t_begin'] = t2t_begin
                    else:
                        df['t2t_begin'] = np.nan
                    # t2t end
                    if 't2t_end' in df.columns:
                        if line['t2t_end'] == '-':
                            df.loc[idx, 't2t_end'] = t2t_end
                    else:
                        df['t2t_end'] = np.nan
                    # hg19 chrom
                    if 'hg19_chrom' in df.columns:
                        if line['hg19_chrom'] == '-':
                            df.loc[idx, 'hg19_chrom'] = hg19_chr
                    else:
                        df['hg19_chrom'] = np.nan
                    # hg19 transcript
                    if 'hg19_transcript' in df.columns:
                        if line['hg19_transcript'] == '-':
                            df.loc[idx, 'hg19_transcript'] = hg19_transcript
                    else:
                        df['hg19_transcript'] = np.nan
                    # hg19 begin
                    if 'hg19_begin' in df.columns:
                        if line['hg19_begin'] == '-':
                            df.loc[idx, 'hg19_begin'] = hg19_begin
                    else:
                        df['hg19_begin'] = np.nan
                    # hg19 end
                    if 'hg19_end' in df.columns:
                        if line['hg19_end'] == '-':
                            df.loc[idx, 'hg19_end'] = hg19_end
                    else:
                        df['hg19_end'] = np.nan

                df['gene_name'] = gene_name
                df['ncbi_id'] = ncbi_id
                df['ncbi_url'] = url
                # 调整顺序
                col_list = ['gene_name', 'ncbi_id',
                            'hg38_chrom', 'hg38_transcript', 'hg38_begin', 'hg38_end',
                            't2t_chrom', 't2t_transcript', 't2t_begin', 't2t_end',
                            'hg19_chrom', 'hg19_transcript', 'hg19_begin', 'hg19_end', 'ncbi_url']
                dataframe = df[col_list]
                return dataframe
        else:
            return df
