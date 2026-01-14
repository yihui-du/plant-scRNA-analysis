# Git仓库上传完整指南

## 重要前置准备

### 步骤1：安装Git（如果未安装）

#### Windows用户
1. 访问 https://git-scm.com/download/win
2. 下载并运行安装程序
3. 使用默认设置完成安装
4. 重启PowerShell或命令行

#### macOS用户
```bash
brew install git
```

#### Linux用户
```bash
sudo apt-get install git  # Ubuntu/Debian
sudo yum install git       # CentOS/RHEL
```

### 步骤2：配置Git

```bash
# 设置用户名
git config --global user.name "Your Name"

# 设置邮箱
git config --global user.email "your.email@example.com"

# 验证配置
git config --list
```

### 步骤3：创建GitHub账户（如果没有）
1. 访问 https://github.com/signup
2. 填写用户名、邮箱和密码
3. 验证邮箱
4. 创建账户

---

## 创建GitHub仓库

### 步骤1：在GitHub上创建新仓库
1. 登录GitHub后，点击右上角"+"图标
2. 选择"New repository"
3. 填写以下信息：
   - Repository name: `plant-scRNA-analysis`
   - Description: `Plant single-cell RNA-seq downstream analysis pipeline`
   - Visibility: 选择 "Public"（如果想分享）或 "Private"（仅自己访问）
4. **不要勾选** "Initialize this repository with:"
5. 点击"Create repository"

### 步骤2：获取远程仓库链接
创建仓库后，GitHub会显示设置页面，找到如下形式的链接：
```
https://github.com/YOUR-USERNAME/plant-scRNA-analysis.git
```
复制这个链接

---

## 本地上传流程

### 步骤1：初始化本地仓库

进入项目目录：
```bash
cd d:\test\plant-scRNA-analysis
```

初始化Git仓库：
```bash
git init
```

### 步骤2：添加文件到暂存区
```bash
# 添加所有文件
git add .

# 或添加特定文件
git add README.md
git add scripts/
```

### 步骤3：创建初始提交
```bash
git commit -m "Initial commit: Project structure and documentation"
```

### 步骤4：添加远程仓库
将 `YOUR-USERNAME` 替换为你的GitHub用户名：
```bash
git remote add origin https://github.com/YOUR-USERNAME/plant-scRNA-analysis.git
```

### 步骤5：推送到GitHub
```bash
# 如果是第一次，可能需要指定分支
git branch -M main
git push -u origin main

# 之后可以直接使用
git push
```

### 步骤6：验证
1. 打开浏览器访问你的仓库: https://github.com/YOUR-USERNAME/plant-scRNA-analysis
2. 检查所有文件是否正确上传

---

## 后续维护

### 上传新代码

当你有新的代码文件要上传时：

```bash
# 1. 复制代码文件到相应目录
# 2. 查看改动
git status

# 3. 添加改动
git add scripts/01_fastq_conversion/your_new_script.py

# 4. 提交
git commit -m "Add FASTQ conversion script"

# 5. 推送到GitHub
git push
```

### 记录问题和解决方案

编辑 `issues/CHANGELOG.md` 和 `docs/troubleshooting.md`：

```bash
git add issues/CHANGELOG.md
git commit -m "Update: Record issue #1 and solution"
git push
```

### 常用命令参考

```bash
# 查看状态
git status

# 查看修改内容
git diff

# 查看提交历史
git log --oneline

# 查看远程仓库信息
git remote -v

# 撤销未提交的改动
git checkout -- filename

# 查看某个文件的历史
git log -- filename

# 创建新分支（用于新功能开发）
git checkout -b feature/new-feature

# 切换分支
git checkout main

# 合并分支
git merge feature/new-feature

# 删除分支
git branch -d feature/new-feature
```

---

## 使用SSH密钥（更安全）

### 生成SSH密钥

```bash
ssh-keygen -t rsa -b 4096 -C "your.email@example.com"
```

按提示操作（可直接按Enter使用默认设置）

### 添加SSH密钥到GitHub

1. 复制公钥内容：
   ```bash
   cat ~/.ssh/id_rsa.pub  # macOS/Linux
   # 或在Windows中用记事本打开 C:\Users\YourUsername\.ssh\id_rsa.pub
   ```

2. 访问 https://github.com/settings/keys
3. 点击"New SSH key"
4. 粘贴公钥内容
5. 点击"Add SSH key"

### 使用SSH克隆/推送

```bash
# 使用SSH而不是HTTPS
git remote add origin git@github.com:YOUR-USERNAME/plant-scRNA-analysis.git

# 或修改现有的远程仓库
git remote set-url origin git@github.com:YOUR-USERNAME/plant-scRNA-analysis.git
```

---

## GitHub Actions自动化（可选）

在项目根目录创建 `.github/workflows/` 目录，可以自动化：
- 代码检查
- 测试运行
- 文档生成

```yaml
# .github/workflows/tests.yml
name: Run Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - name: Set up Python
        uses: actions/setup-python@v2
        with:
          python-version: 3.9
      - name: Install dependencies
        run: pip install -r requirements.txt
      - name: Run tests
        run: python -m pytest
```

---

## 故障排除

### 问题1：git命令未找到
**解决**：检查Git是否正确安装
```bash
git --version
```

### 问题2：Permission denied
**解决**：检查SSH权限
```bash
chmod 600 ~/.ssh/id_rsa
chmod 644 ~/.ssh/id_rsa.pub
```

### 问题3：认证失败
**解决**：
- 检查用户名和密码
- 或使用个人访问令牌（Personal Access Token）代替密码

### 问题4：大文件无法推送
**解决**：
- 使用Git LFS (Large File Storage)
- 或将大文件存储在data目录并在.gitignore中指定

```bash
# 安装Git LFS
git lfs install

# 追踪大文件
git lfs track "*.h5ad"
git add .gitattributes
git commit -m "Add LFS tracking"
```

---

## 最佳实践

1. **经常提交**：定期保存工作
2. **写好提交信息**：描述改动的内容
3. **使用分支**：为新功能创建分支
4. **定期拉取**：保持本地版本最新
5. **备份重要数据**：大数据文件不要放在Git中

---

## 推荐资源

- [Git官方文档](https://git-scm.com/doc)
- [GitHub帮助中心](https://docs.github.com)
- [Pro Git书籍](https://git-scm.com/book/en/v2)（免费在线版）
- [Interactive Git Tutorial](https://learngitbranching.js.org/)

---

**最后更新**: 2026年1月14日
