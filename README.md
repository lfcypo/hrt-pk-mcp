# HRT PK MCP

激素替代治疗 (HRT) 药代动力学血药浓度监测 MCP 工具

Agent 可通过工具记录给药时间，并随时查询根据 PK 算法拟合的血药浓度

参考 [HRT-Recorder-PKcomponent-Test](https://github.com/LaoZhong-Mihari/HRT-Recorder-PKcomponent-Test) 算法与数据拟合实现

支持雌二醇 (Estradiol) 和睾酮 (Testosterone) 的多种给药途径

## 快速开始

```bash
git clone https://github.com/lfcypo/hrt-pk-mcp
uv sync

uv run hrt-pk-mcp
# 或
python -m pk_mcp
```

## MCP 工具

| 工具                    | 说明             |
|-----------------------|----------------|
| `log_dose`            | 记录一次给药事件       |
| `query_concentration` | 查询指定时间点的预测血药浓度 |
| `list_dose_events`    | 列出所有已记录的给药事件   |
| `remove_event`        | 按 ID 删除指定事件    |
| `clear_all_events`    | 清空所有给药事件       |

### log_dose 记录一次给药事件

记录一次激素用药事件。每个事件会被持久化存储，并在后续所有浓度查询中纳入计算。

| 参数                        | 类型     | 说明                                                                           |
|---------------------------|--------|------------------------------------------------------------------------------|
| `compound`                | str    | 药物化合物：`E2`, `EB`, `EV`, `EC`, `EN`, `T`, `TC`, `TE`, `TU`                    |
| `route`                   | str    | 给药途径：`injection`, `gel`, `patch_apply`, `patch_remove`, `oral`, `sublingual` |
| `time_h`                  | float  | 给药时间 (距参考起点的小时数)                                                             |
| `dose_mg`                 | float  | 剂量（活性激素当量, 毫克）                                                               |
| `release_rate_ug_per_day` | float? | 贴片零阶释放速率（µg/d），仅 `patch_apply` 使用                                            |
| `area_cm2`                | float? | 凝胶涂抹面积（cm²），仅 `gel` 使用                                                       |
| `sublingual_theta`        | float? | 舌下含服中经黏膜快通路吸收比例 [0-1]，仅 `sublingual` 使用                                      |

### query_concentration 查询指定时间点的预测血药浓度

| 参数               | 类型    | 说明                                |
|------------------|-------|-----------------------------------|
| `time_h`         | float | 查询时间点（距参考起点的小时数）                  |
| `hormone`        | str   | 目标激素：`estradiol` 或 `testosterone` |
| `body_weight_kg` | float | 体重（kg，默认 70.0；影响分布体积 Vd）          |

## 数据存储

给药事件以 JSON 格式持久化存储在项目目录下的 `.hrt-pk-data/dose_events.json` 文件中

可跨会话检查、备份或迁移

## 药代动力学算法

本项目药代动力学算法与常数参考自 [HRT-Recorder-PKcomponent-Test](https://github.com/LaoZhong-Mihari/HRT-Recorder-PKcomponent-Test)

关于详细的PK算法、完整参数推导、校准来源等，请参阅 [HRT-Recorder](https://github.com/LaoZhong-Mihari/HRT-Recorder-PKcomponent-Test/tree/main/HRT-Recorder)
与 [pk_research](https://github.com/LaoZhong-Mihari/HRT-Recorder-PKcomponent-Test/tree/main/pk_research)

### 支持的化合物与给药途径

#### 雌二醇 (Estradiol)

| 化合物      | 注射 | 口服 | 舌下 | 凝胶 | 贴片 |
|----------|----|----|----|----|----|
| E2（雌二醇）  | —  | ✓  | ✓  | ✓  | ✓  |
| EB（苯甲酸酯） | ✓  | —  | —  | —  | —  |
| EV（戊酸酯）  | ✓  | ✓  | ✓  | —  | —  |
| EC（环丙酸酯） | ✓  | —  | —  | —  | —  |
| EN（庚酸酯）  | ✓  | —  | —  | —  | —  |

#### 睾酮 (Testosterone)

| 化合物      | 注射 | 口服 | 凝胶 | 贴片 |
|----------|----|----|----|----|
| T（睾酮）    | —  | —  | ✓  | ✓  |
| TC（环丙酸酯） | ✓  | —  | —  | —  |
| TE（庚酸酯）  | ✓  | —  | —  | —  |
| TU（十一酸酯） | ✓  | ✓  | —  | —  |

### 输出单位

- 雌二醇 pg/mL pmol/L
- 睾酮 ng/dL ng/mL nmol/L

### 迁移的交叉验证

本项目已进行模拟数据交叉验证实验，本项目运行结果与源项目运行结果数值上一致

```bash
python ./test/cross_validate.py
```

## 许可

本项目遵守 MIT License

## 致谢

* PK算法及常数、算法应用 [LaoZhong-Mihari/HRT-Recorder-PKcomponent-Test](https://github.com/LaoZhong-Mihari/HRT-Recorder-PKcomponent-Test)
* PK算法及常数 [TransmtfTeam/Transmtf-HRT-Tracker](https://github.com/TransmtfTeam/Transmtf-HRT-Tracker)
* 算法应用 [mkx173/Featherline](https://github.com/mkx173/Featherline)
* 算法应用 [xunxunProjects/Oyama-s-HRT-Tracker](https://github.com/xunxunProjects/Oyama-s-HRT-Tracker)

## 修正联系

若对项目中的药代动力学算法部分有任何问题或建议，请提出 issue 或联系学术通讯邮箱 lfcypo@gmail.com

## 严重警告

本项目的以教育目的开发，仅建议用于药代动力学模型研究、算法验证、数据分析方面

本项目输出的血药浓度结果来自数学模型拟合与理论参数计算，不等同于真实个体的血药浓度检测结果，不能替代临床实验室检测、药代动力学研究或医学专业评估

**不建议**将本项目用于任何临床医学用途，包括但不限于：

- 指导使用者开始、停止、调整或改变激素治疗方案；
- 决定雌二醇、睾酮或其他性激素药物的剂量；
- 判断药物治疗是否安全、有效或达到治疗目标；
- 替代医生进行处方、诊疗或医学决策；
