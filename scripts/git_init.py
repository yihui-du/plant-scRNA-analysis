#!/usr/bin/env python
"""
Git初始化和上传脚本
用于快速初始化仓库并上传到GitHub
"""

import os
import subprocess
import sys
from pathlib import Path


def run_command(cmd, description):
    """运行命令并显示结果"""
    print(f"\n{'='*60}")
    print(f"📝 {description}")
    print(f"{'='*60}")
    print(f"命令: {cmd}\n")
    
    result = subprocess.run(cmd, shell=True)
    
    if result.returncode != 0:
        print(f"❌ 命令执行失败！")
        return False
    
    print(f"✅ 成功！")
    return True


def init_git_repo():
    """初始化Git仓库"""
    print("\n🚀 开始初始化Git仓库...\n")
    
    # 检查是否已初始化
    if Path('.git').exists():
        print("⚠️  Git仓库已初始化")
        return True
    
    # 初始化仓库
    if not run_command("git init", "初始化Git仓库"):
        return False
    
    # 添加所有文件
    if not run_command("git add .", "添加所有文件"):
        return False
    
    # 首次提交
    if not run_command(
        'git commit -m "Initial commit: Create project structure and documentation"',
        "创建初始提交"
    ):
        return False
    
    return True


def create_readme_for_github():
    """创建GitHub上传说明"""
    content = """# GitHub上传说明

## 在GitHub上创建新仓库后

1. **在GitHub上创建新的仓库** (不要初始化README.md)
   - 访问 https://github.com/new
   - 输入仓库名: `plant-scRNA-analysis`
   - 选择 "Public" 
   - 不要初始化任何文件
   - 点击 "Create repository"

2. **添加远程仓库**
   ```bash
   git remote add origin https://github.com/YOUR-USERNAME/plant-scRNA-analysis.git
   ```

3. **推送到GitHub**
   ```bash
   git branch -M main
   git push -u origin main
   ```

4. **验证**
   - 访问你的GitHub仓库链接检查文件

## 后续更新

编辑文件后：
```bash
git add .
git commit -m "描述你的改动"
git push
```

## 常见命令

```bash
# 查看状态
git status

# 查看提交历史
git log --oneline

# 查看远程仓库
git remote -v

# 更新本地仓库
git pull origin main
```

## 贡献建议

- 在编辑前创建分支: `git checkout -b feature/new-feature`
- 编辑完成后推送: `git push origin feature/new-feature`
- 在GitHub上创建Pull Request

---

**创建时间**: 2026-01-14
"""
    
    with open("GITHUB_SETUP.md", "w", encoding="utf-8") as f:
        f.write(content)
    
    print("✅ 已创建 GITHUB_SETUP.md")


def main():
    """主函数"""
    print("""
    ╔═══════════════════════════════════════════════════════════╗
    ║   植物单细胞测序分析 - Git仓库初始化脚本                    ║
    ║   Plant scRNA-seq Analysis - Git Repository Setup        ║
    ╚═══════════════════════════════════════════════════════════╝
    """)
    
    # 初始化Git仓库
    if not init_git_repo():
        print("\n❌ Git仓库初始化失败！")
        sys.exit(1)
    
    # 创建GitHub设置说明
    create_readme_for_github()
    
    # 显示后续步骤
    print("""
    ╔═══════════════════════════════════════════════════════════╗
    ║   ✅ Git初始化完成！                                        ║
    ╠═══════════════════════════════════════════════════════════╣
    ║   后续步骤：                                               ║
    ║   1. 在 GitHub 上创建新仓库                               ║
    ║   2. 运行以下命令添加远程仓库：                            ║
    ║      git remote add origin https://github.com/YOUR-USERNAME/plant-scRNA-analysis.git
    ║   3. 推送到 GitHub:                                        ║
    ║      git branch -M main                                   ║
    ║      git push -u origin main                              ║
    ║   4. 查看 GITHUB_SETUP.md 了解更多细节                     ║
    ╚═══════════════════════════════════════════════════════════╝
    """)


if __name__ == "__main__":
    main()
