# 双端 FASTQ 质控 API

FastAPI 实现的成对 FASTQ 质量过滤服务：上传 `R1.fastq.gz` / `R2.fastq.gz`
与 Phred+33 质量阈值，后台流式处理，支持查询状态、停止任务、下载过滤后的
成对 gzip 结果。成功任务还可提交**去标识导出**，生成外协分析可用的成对 gzip，
同时在服务端数据目录留下读段级映射表。导出入口同样在 `/docs` 中可见。

## 快速启动

```bash
docker compose up --build
```

启动后：

* Swagger 文档：<http://localhost:8000/docs>
* OpenAPI JSON：<http://localhost:8000/openapi.json>
* 健康检查：<http://localhost:8000/health>

数据目录 `/data`（Compose 卷 `qc-data`）保存上传文件、结果，以及
`/data/exports/<export_id>/` 下的导出产物和映射表。容器重启后内存中的任务/
导出记录会失效；未完成的导出目录会被自动清理，已完成的导出目录（含映射表）
保留。

不装 Docker 时也可本地运行：

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
QC_DATA_DIR=./data uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 接口一览

| 方法 | 路径 | 说明 |
|------|------|------|
| `POST` | `/jobs` | 上传 `r1`、`r2`（multipart `File`）和 `min_quality`，返回 `202` 与任务号 |
| `GET` | `/jobs/{job_id}` | 查询状态、进度、保留/丢弃对数、错误信息、下载链接 |
| `GET` | `/jobs` | 列出所有任务（新到旧） |
| `POST` | `/jobs/{job_id}/stop` | 请求停止；已完成/失败的任务返回 `409` |
| `GET` | `/jobs/{job_id}/download/r1` | 下载过滤后的 R1（仅成功任务） |
| `GET` | `/jobs/{job_id}/download/r2` | 下载过滤后的 R2（仅成功任务） |
| `POST` | `/jobs/{job_id}/exports` | 为成功任务提交去标识导出，返回 `202` 与导出号 |
| `POST` | `/exports` | 同上，任务号放在 JSON 请求体 `{"job_id": "..."}` 中 |
| `GET` | `/jobs/{job_id}/exports` | 列出某任务的导出（新到旧） |
| `GET` | `/exports` | 列出所有导出（新到旧） |
| `GET` | `/exports/{export_id}` | 查询导出状态、进度、导出读段对数和下载链接 |
| `POST` | `/exports/{export_id}/stop` | 请求停止导出；已结束的导出返回 `409` |
| `GET` | `/exports/{export_id}/download/r1` | 下载去标识 R1（仅成功导出） |
| `GET` | `/exports/{export_id}/download/r2` | 下载去标识 R2（仅成功导出） |

`min_quality` 为 0..93 的整数（Phred+33 可表示的范围）。上传字段名固定为
`r1`、`r2`。

### 示例

```bash
# 提交任务
curl -X POST http://localhost:8000/jobs \
  -F 'r1=@samples/mixed_R1.fastq.gz' \
  -F 'r2=@samples/mixed_R2.fastq.gz' \
  -F 'min_quality=20'

# 查询
curl http://localhost:8000/jobs/<job_id>

# 停止
curl -X POST http://localhost:8000/jobs/<job_id>/stop

# 成功后可下载
curl -OJ http://localhost:8000/jobs/<job_id>/download/r1
curl -OJ http://localhost:8000/jobs/<job_id>/download/r2

# 成功任务提交去标识导出
curl -X POST http://localhost:8000/jobs/<job_id>/exports
# 或者
curl -X POST http://localhost:8000/exports \
  -H 'Content-Type: application/json' \
  -d '{"job_id": "<job_id>"}'

# 查询 / 停止导出
curl http://localhost:8000/exports/<export_id>
curl -X POST http://localhost:8000/exports/<export_id>/stop

# 导出成功后下载成对 gzip（文件名 R1.deidentified.fastq.gz / R2...）
curl -OJ http://localhost:8000/exports/<export_id>/download/r1
curl -OJ http://localhost:8000/exports/<export_id>/download/r2
```

### 状态字段

`status` 取值：`pending` → `running` → `succeeded` / `failed`，或
`running`/`pending` → `stopping` → `stopped`。

成功时：

```json
{
  "status": "succeeded",
  "total_pairs": 7,
  "kept_pairs": 4,
  "discarded_pairs": 3,
  "downloads": {
    "r1": "/jobs/<id>/download/r1",
    "r2": "/jobs/<id>/download/r2"
  },
  "error": null
}
```

导出成功时（`pairs_exported` 与 `total_pairs` 同为导出读段对数）：

```json
{
  "export_id": "<export_id>",
  "job_id": "<source_job_id>",
  "source_job_id": "<source_job_id>",
  "status": "succeeded",
  "pairs_processed": 4,
  "pairs_exported": 4,
  "total_pairs": 4,
  "downloads": {
    "r1": "/exports/<export_id>/download/r1",
    "r2": "/exports/<export_id>/download/r2"
  },
  "error": null
}
```

失败时 `kept_pairs`/`discarded_pairs`（或导出的 `pairs_exported`）/`downloads`
为 `null`，`error.code` 为下列之一：

| code | 含义 |
|------|------|
| `pair_name_mismatch` | 去掉 `/1`、`/2` 后缀后两端名称不同 |
| `invalid_name_suffix` | 名称未以 `/1`（R1）或 `/2`（R2）结尾 |
| `pair_count_mismatch` | 两端记录条数不同 |
| `seq_qual_length_mismatch` | 序列与质量行长度不同 |
| `malformed_record` | 记录不是完整四行、缺少 `@`/`+`、行超长等 |
| `invalid_quality_char` | 质量字符超出 Phred+33 可打印范围（ASCII 33..126） |
| `gzip_corrupt` | gzip 截断、CRC 错误、尾部垃圾、非 gzip 文件或空文件 |
| `mapping_write_failed` | （导出）映射表写入或落盘失败 |
| `write_failed` | （导出）成对 gzip 输出写入或收尾失败 |
| `publish_failed` | （导出）临时文件改名发布失败 |
| `internal_error` | 未预期的内部错误 |

### 去标识导出规则

导出只接受**已成功**的质控任务，后台逐对读取该任务已发布的
`R1.kept.fastq.gz` / `R2.kept.fastq.gz`，因此导出的就是质控保留的读段对，
顺序不变。

1. **新编号**：每一读段对获得同一个新编号。导出的头部为
   `@DX<导出号前 12 位>_<12 位十进制序号>/1`（R1）和 `.../2`（R2），
   两端编号部分完全相同；新编号只由导出号与对序号决定，不含任何原名称信息。
2. **去除旧标识**：原头部整行（含空格后的注释）被替换，`+` 行附带的名称
   也被替换为裸 `+`；导出文件的标识字段不会保留旧名称。
3. **序列与质量**：逐字节原样复制（仅行尾统一为 `\n`）。
4. **映射表**：逐对写入 `/data/exports/<export_id>/mapping.jsonl`，每行一个
   JSON 对象：
   `{"pair_number": 1, "new_id": "DX..._000000000001", "r1_header": "@原R1完整头部", "r2_header": "@原R2完整头部"}`
   （`r1_header`/`r2_header` 为去掉行尾后的原始完整头部）。映射表**只存于
   服务端数据目录**，没有任何 API 提供读取或下载；导出接口的响应与下载文件
   均不含映射内容。
5. **成对结果**：成功后下载 R1、R2 两个 gzip，响应中的
   `pairs_exported`（亦作 `total_pairs`）为读段对数。
6. **流式处理**：逐对读取、逐对写出，不整份载入内存；每 256 对更新进度。
7. **失败清理**：源文件损坏（gzip 截断/CRC 错误等）、映射表写入失败、任一
   输出 gzip 写入/收尾失败、发布改名失败，均会删除本次导出的全部产物
   （输出、映射表、`.part` 文件、整个导出目录）后才报告 `failed`；源任务的
   质控结果不受影响。
8. **停止**：导出可查询、可停止；一旦接受停止即进入 `stopping`，等到后台
   写入退出并清理完整个导出目录后才变为 `stopped`，之后**不可能**再转为
   `succeeded`。导出运行期间的其它任务（质控任务）查询与停止照常响应。


## 处理规则

1. 每条读段恰好四行：`@名称`、序列、`+`、质量。
2. R1 名称必须以 `/1` 结尾，R2 必须以 `/2` 结尾；去掉后缀后逐对核对，
   必须完全一致。名称后的注释（空格分隔）会被忽略。
3. 序列与质量长度必须一致。
4. 质量按 Phred+33 解释；只要任一端出现低于阈值的碱基，整对丢弃。
5. 保留的读段按原顺序写出，输出名称沿用原始四行（仅补一个 `\n`）。
6. gzip 损坏、条数不符、名称不符、长度不符均导致任务 `failed`，不产生任何
   可下载文件。
7. 流式处理：逐对读取、逐对写出，不会整份载入内存；每 256 对更新一次进度。
   处理期间状态查询与停止请求由独立线程和事件立即响应。
8. 停止与完成互斥：一旦接受停止，任务不可能转为成功；完成后再次停止返回
   `409`。`stopping` 只在工作线程退出且临时文件（上传文件、`.part` 中间
   输出）清空后，才会变为 `stopped`。
9. 只有 `succeeded` 任务提供下载；失败或停止均返回 `409`，不会给出半成品。
   输出先在 `.part` 文件中完整写入，成功后才原子改名发布。

## 本地样例

`tools/make_samples.py` 可重新生成 `samples/` 下的样例：

| 文件 | 内容 |
|------|------|
| `good_R1/R2.fastq.gz` | 3 对全部通过的读段 |
| `bad_R1/R2.fastq.gz` | 3 对全部因低质量被丢弃 |
| `mixed_R1/R2.fastq.gz` | 7 对，前 3 与最后 1 对保留，中间 3 对丢弃（阈值 20） |
| `truncated_R1.fastq.gz` + `truncated_R2.fastq.gz` | R2 的 gzip 流末尾被截断，任务必定失败为 `gzip_corrupt` |

```bash
python tools/make_samples.py
```

## 测试

```bash
pip install -r requirements.txt pytest httpx
python -m pytest tests -q
```

测试覆盖过滤/顺序/阈值边界、各类强制失败、停止语义（含停止抢先于完成/
失败、清理时序）、处理期间查询响应性、流式读取行数上限，以及去标识导出的
编号替换、映射表内容与服务端隔离、各类写入/发布失败清理和导出停止语义。

## 说明与限制

* 单进程内存任务表，符合“单个 API 进程管理任务”的范围；不要用多 worker
  启动（`--workers 1`），否则任务/导出查询会落到不同进程。Compose 默认
  单进程。
* 每个上传文件默认上限 20 GiB，可用 `QC_MAX_UPLOAD_BYTES` 调整。
* 后台线程为守护线程；进程关闭时会给所有任务/导出发停止信号并等待清理。
* 导出映射表随导出目录保存在数据卷 `/data/exports/<export_id>/mapping.jsonl`，
  需要长期追溯时请备份该卷；任务或导出在内存中的记录会在进程重启后消失，
  但已完成导出的目录不会被启动清理删除。若需回收磁盘，可删除
  `/data/exports/<export_id>/`（映射表也会一并删除）。

