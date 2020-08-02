# %%
import re
import string
from itertools import compress
from pypinyin import lazy_pinyin


# %%
def _style_a_name(name_with_space: str):
    """
    对单个英文名字进行格式化。空格删掉，family name 后加','，后面缩写为一个字母的first name 后加'.'

    parameter:
    ----------
        name_with_space: str
            带空格的单个名字

    return:
        str。格式化后的单个名字
    """

    name_parts = name_with_space.split(' ')
    # 如果只有一项，就直接返回
    if len(name_parts) == 1:
        return name_with_space
    replaced_parts = [
        part + ',' if len(part) > 1 else part + '.' for part in name_parts
    ]
    return ''.join(replaced_parts)


# %%
def format_a_line(line: str):
    """
    对于一行的reference 进行格式化。如果元素不足导致失败将直接返回原reference 值。

    Parameter:
    ----------
        line: str
            一行自动生成的国标refernce 行

    Return:
    -------
        str 整理过格式的同行reference
    """
    # "(.+)，\. *(.+?)\. *(.+?), *(\d{4})(.+)"
    # (.+)，\. *(.+?)\. *(.+)
    mat = re.search(r"(.+)，\. *(.+?)\. *(.+?), *(\d{4}),*(.+)", line)
    try:
        names, title, magzine, year, others = mat.groups()
    except AttributeError:
        # 如果正则解析失败，则直接返回原数据
        return line

    # 把名字拆成list
    names_list = names.title().split(', ')[:3]
    styled_names = [_style_a_name(name) for name in names_list]

    # 对中文与英文分开处理
    if names[0] in string.ascii_uppercase:
        # 合并名字字符串，并在最后加上', et al'
        final_names_str = ", ".join(styled_names)
        if len(styled_names) > 3:
            final_names_str += ', et al'

        # 英文格式的最终板式
        ref_result = '{name}, {year}, "{title}", {magzine}, {others}'.format(
            name=final_names_str,
            year=year,
            title=title.title(),
            magzine=magzine,
            others=others.strip())
    else:
        final_names_str = "、".join(styled_names)
        if len(styled_names) > 3:
            final_names_str += '，等'

        # 中文格式的最终结果
        ref_result = '{name}. {title}. {magzine}，{year}, {others}'.format(
            name=final_names_str,
            title=title,
            magzine=magzine,
            year=year,
            others=others.strip())
    return ref_result


# %%
if __name__ == "__main__":
    with open("text/gb_refer.txt", "r") as f:
        raw_txt = f.read()
    raw_lines = raw_txt.splitlines()

    # 对每一行的ref 进行格式化
    formatted_refs = [format_a_line(raw_ref) for raw_ref in raw_lines]

    # 中英文分开排序并合并
    en_mask = [ref[0] in string.ascii_letters for ref in formatted_refs]
    zh_mask = [not en for en in en_mask]
    zh_refs = list(compress(formatted_refs, zh_mask))
    en_refs = list(compress(formatted_refs, en_mask))

    zh_sorted = sorted(zh_refs, key=lambda ref: lazy_pinyin(ref[0]))
    en_sorted = sorted(en_refs)
    all_sorted = zh_sorted + en_sorted

    # 加序号
    serial_num = ['[{}] '.format(num) for num in range(1, len(all_sorted) + 1)]
    all_with_serial = [num + ref for num, ref in zip(serial_num, all_sorted)]

    # 转换为本文
    formatted_refs_txt = '\n'.join(all_with_serial)
    result_refs_txt = re.sub(r', *\.', '.', formatted_refs_txt)

    with open("text/formatted_refs.txt", 'w') as f:
        f.write(result_refs_txt)
