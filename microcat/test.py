#!/usr/bin/env python
import argparse


def add(args):
    result = args.num1 + args.num2
    print(f"{args.num1} + {args.num2} = {result}")

def subtract(args):
    result = args.num1 - args.num2
    print(f"{args.num1} - {args.num2} = {result}")

def multiply(args):
    result = args.num1 * args.num2
    print(f"{args.num1} * {args.num2} = {result}")

def divide(args):
    result = args.num1 / args.num2
    print(f"{args.num1} / {args.num2} = {result}")

def main():
    parser = argparse.ArgumentParser(description="一个支持多种操作的计算器")

    # 创建子命令解析器
    subparsers = parser.add_subparsers(help="支持的命令")

    # 添加 add 子命令
    parser_add = subparsers.add_parser("add", help="执行加法操作")
    parser_add.add_argument("num1", type=int, help="第一个数字")
    parser_add.add_argument("num2", type=int, help="第二个数字")
    parser_add.set_defaults(func=add)

    # 添加 subtract 子命令
        # 添加 subtract 子命令
    parser_subtract = subparsers.add_parser("subtract", help="执行减法操作")
    parser_subtract.add_argument("num1", type=int, help="第一个数字")
    parser_subtract.add_argument("num2", type=int, help="第二个数字")
    parser_subtract.set_defaults(func=subtract)

    # 添加 multiply 子命令
    parser_multiply = subparsers.add_parser("multiply", help="执行乘法操作")
    parser_multiply.add_argument("num1", type=int, help="第一个数字")
    parser_multiply.add_argument("num2", type=int, help="第二个数字")
    parser_multiply.set_defaults(func=multiply)

    # 添加 divide 子命令
    parser_divide = subparsers.add_parser("divide", help="执行除法操作")
    parser_divide.add_argument("num1", type=int, help="第一个数字")
    parser_divide.add_argument("num2", type=int, help="第二个数字")
    parser_divide.set_defaults(func=divide)

    # 解析参数并执行相应的函数
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
