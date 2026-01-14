# GitHub上传说明

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
