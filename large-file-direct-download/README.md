# large-file-direct-download (`direct-download`)

在 **不改动 Claude Code / shell 代理** 的前提下，让指定的大文件下载任务走本地网络 DIRECT，
不消耗 Clash / Mihomo 代理节点流量。单文件 Python 3 标准库脚本，依赖 `ss`、`ip`、`curl`，`fetch` 需要 `aria2c`。

## 原理（为什么不改 Mihomo 规则）

本机实测（2026-11-10）：FlClash (Mihomo) `mode: global`、`tun.enable: false`、无 tun 网卡、无额外 `ip rule`，
默认路由直接走 `eno1`。此时代理只通过 `http(s)_proxy` 等环境变量 / 应用自身设置生效，
**不连 127.0.0.1:7890 的进程本来就直接出网**。所以最小方案是：

1. **只给下载子进程一个干净环境**：删除 `http/https/ftp/all/no_proxy`（大小写）、设 `NO_PROXY='*'`、
   `WGETRC` 指向 `use_proxy = off`、`GIT_CONFIG_*` 令 `http.proxy=''`；`fetch` 还给 aria2c 加 `--no-conf --all-proxy=`。
   调用者（Claude Code）的环境不变。子进程派生的所有进程（organise-sra → kingfisher → aria2c/prefetch…）继承干净环境。
2. **socket 级证明**：运行期间每 0.5 s 用 `ss -tinp` 扫描整个进程树（同 session + 所有后代）的 TCP 连接：
   - 连到本地代理端口 → `PROXY` 泄漏 → 默认立即终止任务（退出码 3）；
   - 远端地址用 `ip route get` 取出口网卡，tun/wg 类网卡 → `TUN` 泄漏；
   - 其余记为 `DIRECT`，并统计每个连接收到的字节数；
   - 同时统计代理端口上所有程序的字节增量，作为"代理流量计数"参考（external-controller 关闭，无法读 Mihomo API）。
3. **TUN 开启时拒绝执行**（检测到 tun 类网卡、Mihomo `tun.enable: true`、非默认 `ip rule` 或默认路由走 tun），
   因为那时剥离环境变量不再足够，工具不会假装成功。

> "DIRECT" 的含义：**未经过 Mihomo 代理节点，且从本机物理网卡（eno1）出去**。上游路由器 / 网关层面的 VPN 本机无法看见。
> 隐私：不访问任何 IP 查询服务，日志不记录本机地址，只记录远端服务器地址。

## 用法

```bash
direct-download check                      # 本地环境审计（代理变量、各工具 rc、Mihomo 配置、TUN/路由）
direct-download selftest                   # 真实小文件测试: curl/wget/aria2c/python 直连 + 经代理的反向对照
direct-download size URL...                # 直连 HEAD：大小、是否支持断点续传

# 任意命令（SRA / organise-sra / kingfisher / prefetch / wget / python 脚本 …）
direct-download run -- organise-sra download ...
direct-download run -- prefetch SRR12118866 --max-size 50G
direct-download run --log-dir out/ -- python my_download.py

# 托管下载（aria2c：断点续传、重试、并行、限速、校验、跳过已完成、报告）
direct-download fetch -o refs/ URL1 URL2 --md5 <md5>      # --md5/--sha256 仅适用单个 URL
direct-download fetch -o refs/ -i manifest.tsv -j 3 -x 4 --limit 30M
```

`manifest.tsv`（表头可选）：`url  [out]  [md5|sha256]  [size]`，例如

```
url	md5
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../x_genomic.fna.gz	c13d459b5caa702ff7e1f26fe44b8ad7
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

### fetch 的完整性规则（不会在没有证据时标记成功）

| 状态 | 条件 |
|---|---|
| `VERIFIED_CHECKSUM` | md5/sha256 一致（aria2c 下载时校验，结束后再算一次） |
| `VERIFIED_SIZE` | 大小 = 服务器 Content-Length / manifest，且压缩文件通过 `gzip/bzip2/xz/zstd -t` |
| `UNVERIFIED` | 无校验值且服务器不报大小 → 退出码 4（`--allow-unverified` 可放行） |
| `FAILED` / `INCOMPLETE` / `MISSING` | 退出码 1；`.aria2` 控制文件保留，重跑即续传 |
| `SKIPPED_*` | 已存在且已验证，不重复下载 |

已存在但无法验证的旧文件会改名为 `*.bad-<时间>`（不删除），然后重新下载。
服务器不支持 Range 时会在日志里说明"重新下载而非续传"。
日志与报告：`<outdir>/.direct-download/{fetch,aria2,report}-<时间>.{log,json,tsv}`。

### 退出码

`0` 成功 · `1` 下载/校验失败 · `2` 环境不能保证直连（TUN 等）· `3` 检测到代理泄漏或未能证明直连 · `4` 有 UNVERIFIED · `5` 主机无法直连（不会回退到代理）

## 限制（如实说明）

- 只监测 **TCP**；UDP/QUIC 下载不在监测范围（wget/curl/aria2c/prefetch/ascp-over-TCP 均为 TCP；Aspera 的 FASP 是 UDP，
  它不走 HTTP 代理，但本工具不能给出 socket 证明）。
- 采样间隔 0.5 s，短于采样间隔的小连接可能看不到；大文件连接是长连接，一定会被看到。`DIRECT=0` 时报告"未证明"而非成功。
- 有些站点直连不可达或很慢；`fetch` 预检直连 TCP 不通就停止（退出码 5），**不会自动改走代理**。
- 如果以后开启 Mihomo TUN：工具会拒绝运行。届时可选方案（都需要你同意）：Mihomo 规则模式下加
  `PROCESS-PATH`/`UID` DIRECT 规则（需 `find-process-mode` ≠ off，或专用下载用户 + `UID,<uid>,DIRECT`），然后用 `--allow-tun`；
  监测器仍会根据出口网卡把走 tun 的连接标为泄漏。

## 安装 / 回滚

```bash
cd Temp_scripts/large-file-direct-download
ln -sf "$PWD/direct-download" ~/.local/bin/direct-download
mkdir -p ~/.claude/skills/large-file-direct-download
cp skill/SKILL.md ~/.claude/skills/large-file-direct-download/      # Claude Code skill (触发条件与调用步骤)
direct-download check && direct-download selftest                  # 安装后先做真实小文件直连测试
```

回滚：`rm ~/.local/bin/direct-download; rm -r ~/.claude/skills/large-file-direct-download`，
并删掉 `~/.claude/CLAUDE.md` 里的 "Large downloads" 一行（如有）。本工具未修改任何系统网络配置、Mihomo/FlClash 配置或代理变量。
