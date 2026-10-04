# LeechKEM 独立输入与完整流程运行时间基准（论文一致版 v4）

本修正版严格采用论文中的
`[g_tau(b_g)]_j = sum_r 2^(r+1) b_{g,24r+j}`。程序除检查编解码往返外，
还会直接检查每个 `g_tau` 输入比特对应的输出坐标。

本版本在精确实现的基础上改进了实验设计：每个方案默认预先生成
31 组相互独立的 `(pk, sk, ct, ss)` 测试数据，合法和非法 Decaps
计时轮换使用这些数据，而不是反复使用同一个密文；此外单独测量
`KeyGen + Encaps + 合法 Decaps` 的完整 KEM 流程。

## 在 VS Code 的 WSL 终端运行

进入解压后的目录，执行：

```bash
chmod +x run_full_benchmark.sh
./run_full_benchmark.sh
```

正式命令默认使用 31 组独立数据、31 个计时批次，并把每批运行时间
校准到至少 20 ms。第一次只想检查安装是否正确，可以运行：

```bash
./run_full_benchmark.sh 5 1
```

短命令生成的数据只用于检查程序，不能作为论文正式结果。正式运行时
建议连接电源、关闭重负载程序，并使用正常或高性能电源模式。

## 五类计时结果

- `keygen`：完整密钥生成。
- `encaps`：完整封装，包括消息和盐生成、哈希、PKE 加密及最终共享密钥派生。
- `decaps_valid`：合法密文的完整解封装，包括 Leech 解码、重加密检查和密钥选择。
- `decaps_invalid`：对应合法密文翻转一位后的完整隐式拒绝路径。
- `full_cycle`：同一次计时调用中依次重新执行 KeyGen、Encaps 和合法 Decaps。

单项 Decaps 计时不包含测试数据的预生成时间。`full_cycle` 则确实把三步
全部放入计时区间，并且每次调用都会生成新的密钥和密文。因此二者回答
不同问题，不应把单项中位数直接相加代替 `full_cycle`。

## 需要发回的结果

正式运行结束后，请发回：

- `results/runtime_results.csv`
- `results/environment.txt`

程序还会自动生成可粘贴到论文中的：

- `results/runtime_results.tex`
- `results/runtime_results.md`
- `results/runtime_ratios.tex`

CSV 中的 `samples` 表示计时批次数，`case_pool_size` 表示独立测试数据组数；
默认情况下二者均为 31。
